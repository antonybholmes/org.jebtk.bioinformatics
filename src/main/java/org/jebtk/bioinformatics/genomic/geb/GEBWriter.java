package org.jebtk.bioinformatics.genomic.geb;

import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.bioinformatics.genomic.TagType;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterHashMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.collections.UniqueArrayListCreator;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.json.Json;
import org.jebtk.core.json.JsonObject;
import org.jebtk.core.sys.SysUtils;
import org.jebtk.core.text.TextUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GEBWriter {
  private final byte[] BUFFER = new byte[GEBReader.MAX_VARCHAR_LENGTH];

  // private char[] mCBuf;

  // private byte[] mN;

  // private byte[] mMask;

  // private byte[] mDna;

  private static final Logger LOG = LoggerFactory.getLogger(GEBWriter.class);

  private int mWindow;

  private boolean mRadixMode;

  private Path mDir;

  private String mPrefix;

  public GEBWriter(Path dir, String prefix, Genome genome, int window) {
    this(dir, prefix, genome, window, true);
  }

  public GEBWriter(Path dir, String prefix, Genome genome, int window, boolean radixMode) {
    mDir = dir;
    mPrefix = prefix;
    mWindow = window;
    mRadixMode = radixMode;
  }

  public <T extends GenomicElement> void write(Collection<T> elements) throws IOException {
    // Path mDir = file.toAbsolutePath().getParent();

    write(GenomicElement.toMap(elements));
  }

  public <T extends GenomicElement> void write(IterMap<Chromosome, Set<T>> elements) throws IOException {
    // Path mDir = file.toAbsolutePath().getParent();

    IterTreeMap<String, Integer> sizeBytesMap = new IterTreeMap<String, Integer>();

    IterTreeMap<Integer, Integer> intSizeBytesMap = new IterTreeMap<Integer, Integer>();

    IterTreeMap<Double, Integer> doubleSizeBytesMap = new IterTreeMap<Double, Integer>();

    IterTreeMap<String, Integer> offsetBytesMap = new IterTreeMap<String, Integer>();

    IterTreeMap<Integer, Integer> intOffsetBytesMap = new IterTreeMap<Integer, Integer>();

    IterTreeMap<Double, Integer> doubleOffsetBytesMap = new IterTreeMap<Double, Integer>();

    IterMap<GenomicElement, Integer> elementOffsetBytes = new IterTreeMap<GenomicElement, Integer>();

    //
    // Determine gene addresses
    //

    int offset = 0;

    for (Entry<Chromosome, Set<T>> item : elements) {
      Chromosome chr = item.getKey();

      Set<T> features = elements.get(chr);

      for (GenomicElement e : features) {
        elementOffsetBytes.put(e, offset);

        offset += elementSizeBytes(e);
      }

      // Assemble all the tags
      tagSizesBytes(features, sizeBytesMap, intSizeBytesMap, doubleSizeBytesMap, offsetBytesMap, intOffsetBytesMap,
          doubleOffsetBytesMap);
    }

    //
    // Now process each chr in turn for the bin files
    //

    for (Entry<Chromosome, Set<T>> item : elements) {
      Chromosome chr = item.getKey();

      Set<T> features = item.getValue();

      DefaultTreeMap<Integer, List<GenomicElement>> binsMap = DefaultTreeMap
          .create(new UniqueArrayListCreator<GenomicElement>());

      int minBin = Integer.MAX_VALUE;
      int maxBin = Integer.MIN_VALUE;

      // Find the max bin in the list
      for (GenomicElement e : features) {

        int bs = e.getStart() / mWindow;
        int be = e.getEnd() / mWindow;

        binsMap.get(bs).add(e);
        binsMap.get(be).add(e);

        minBin = Math.min(minBin, bs);
        maxBin = Math.max(maxBin, be);
      }

      binsMap.setAutoCreate(false);

      // Bins are collections of gene addresses
      IterMap<Integer, Integer> binSizeBytes = new IterTreeMap<Integer, Integer>();

      // int binsWidthBytes = 0;

      for (int b = minBin; b <= maxBin; ++b) {
        // size of bin is number of items + n addresses of items

        int n;

        if (binsMap.containsKey(b)) {
          n = binsMap.get(b).size();
        } else {
          n = 0;
        }

        int s = GEBReader.INT_BYTES * (1 + n);

        binSizeBytes.put(b, s);
      }

      writeBTree(chr, binsMap, binSizeBytes, elementOffsetBytes);

      // writeBins(chr, binsMap, binSizeBytes, elementOffsetBytes);

      // break;
    }

    //
    // Write out the elements
    //

    writeElements(elements, sizeBytesMap, intSizeBytesMap, doubleSizeBytesMap, offsetBytesMap, intOffsetBytesMap,
        doubleOffsetBytesMap);

    // Write the tags file
    writeData(offsetBytesMap, intOffsetBytesMap, doubleOffsetBytesMap);

    if (mRadixMode) {
      writeRadix(elements, elementOffsetBytes);
    }

    // Finally write the index file
    writeIndex();
  }

  private void writeVarchar(String s, DataOutputStream writer) throws IOException {
    int l = Math.min(s.length(), GEBReader.MAX_VARCHAR_LENGTH);
    writer.writeByte(l);
    SysUtils.arraycopy(s.getBytes(StandardCharsets.UTF_8), BUFFER, l);
    writer.write(BUFFER, 0, l);
  }

  /**
   * Return the bytes required to occupy a var char.
   * 
   * @param s
   * @return
   */
  private static int varcharSize(String s) {
    return 1 + Math.min(s.length(), GEBReader.MAX_VARCHAR_LENGTH);
  }

  private <T extends GenomicElement> void writeElements(IterMap<Chromosome, Set<T>> elements,
      IterTreeMap<String, Integer> sizeBytesMap, IterTreeMap<Integer, Integer> intSizeBytesMap,
      IterTreeMap<Double, Integer> doubleSizeBytesMap, IterTreeMap<String, Integer> offsetBytesMap,
      IterTreeMap<Integer, Integer> intOffsetBytesMap, IterTreeMap<Double, Integer> doubleOffsetBytesMap)
      throws IOException {
    int tagsStartBytes = DataReader.HEADER_BYTES_OFFSET + GEBReader.INT_BYTES;

    int intTagsStartBytes = tagsStartBytes + GEBReader.INT_BYTES;

    for (Entry<String, Integer> e : sizeBytesMap) {
      intTagsStartBytes += e.getValue();
    }

    // end of ints plus count for doubles
    int doubleTagsStartBytes = intTagsStartBytes + GEBReader.INT_BYTES;

    for (Entry<Integer, Integer> e : intSizeBytesMap) {
      doubleTagsStartBytes += e.getValue();
    }

    Path file = mDir.resolve(ElementReader.getFileName(mPrefix)); // PathUtils.getPath(mGenome

    DataOutputStream writer = FileUtils.newDataOutputStream(file);

    LOG.info("Writing elements to {}...", file);

    writer.writeInt(GEBReader.CHECK);
    writer.writeByte(GEBReader.VERSION);
    writer.writeInt(mWindow);

    int size = 0;

    for (Entry<Chromosome, Set<T>> item : elements) {
      size += item.getValue().size();
    }

    writer.writeInt(size);

    for (Entry<Chromosome, Set<T>> item : elements) {
      for (GenomicElement e : item.getValue()) {
        // Each gene is given an id to make it easier to check which genes
        // we have found when decoding

        writeElement(writer, e, tagsStartBytes, offsetBytesMap, intTagsStartBytes, intOffsetBytesMap,
            doubleTagsStartBytes, doubleOffsetBytesMap);
      }
    }

    writer.close();
  }

  private static void writeElement(DataOutputStream writer, GenomicElement e, int tagsStartBytes,
      Map<String, Integer> tagOffsetBytes, int intTagsStartBytes, Map<Integer, Integer> intTagOffsetBytes,
      int doubleTagsStartBytes, Map<Double, Integer> doubleTagOffsetBytes) throws IOException {
    // Set all bits to 1
    // writer.write(GEBReader.BLOCK_SEPARATOR);

    writer.writeInt(tagsStartBytes + tagOffsetBytes.get(e.getType().toString()));

    writer.writeInt(tagsStartBytes + tagOffsetBytes.get(e.getChr().toString()));
    // writer.write(GEBReader.getChr(e.mChr));

    writer.writeInt(e.getStart());
    writer.writeInt(e.getEnd());

    // Strand
    writer.write(GEBReader.getStrand(e.getStrand()));

    writer.write(e.getPropertyCount());

    int size = 0;

    for (String name : e.getPropertyNames()) {
      String value = e.getProperty(name);

      writeProperty(writer, name, value, tagsStartBytes, tagOffsetBytes, intTagsStartBytes, intTagOffsetBytes,
          doubleTagsStartBytes, doubleTagOffsetBytes);

      if (++size == GEBReader.MAX_CHILDREN) {
        break;
      }
    }

    writer.write(e.getTagCount());

    size = 0;

    for (String tag : e.getTags()) {
      writeTag(writer, tag, tagsStartBytes, tagOffsetBytes, intTagsStartBytes, intTagOffsetBytes, doubleTagsStartBytes,
          doubleTagOffsetBytes);

      if (++size == GEBReader.MAX_CHILDREN) {
        break;
      }
    }

    // Count children
    size = 0;

    for (Entry<GenomicType, List<GenomicElement>> item : e.getChildren()) {
      size += item.getValue().size();
    }

    writer.writeShort(Math.min(size, GEBReader.MAX_CHILDREN));

    size = 0;

    for (GenomicType type : e.getChildTypes()) {
      for (GenomicElement child : e.getChildren(type)) {
        writeElement(writer, child, tagsStartBytes, tagOffsetBytes, intTagsStartBytes, intTagOffsetBytes,
            doubleTagsStartBytes, doubleTagOffsetBytes);

        if (++size == GEBReader.MAX_CHILDREN) {
          break;
        }
      }
    }
  }

  private static void writeProperty(DataOutputStream writer, String name, Object tag, int tagsStartBytes,
      Map<String, Integer> tagOffsetBytes, int intTagsStartBytes, Map<Integer, Integer> intTagOffsetBytes,
      int doubleTagsStartBytes, Map<Double, Integer> doubleTagOffsetBytes) throws IOException {

    writer.writeInt(tagsStartBytes + tagOffsetBytes.get(name));

    writeTag(writer, tag, tagsStartBytes, tagOffsetBytes, intTagsStartBytes, intTagOffsetBytes, doubleTagsStartBytes,
        doubleTagOffsetBytes);
  }

  private static TagType getTagType(Object tag) {
    TagType t = TagType.TEXT;

    if (tag instanceof Integer) {
      t = TagType.INT;
    } else if (tag instanceof Number) {
      t = TagType.DOUBLE;
    } else {
      t = TagType.TEXT;

      String s = tag.toString();

      if (TextUtils.isInt(s)) {
        t = TagType.INT;
      }

      if (TextUtils.isDouble(s)) {
        t = TagType.DOUBLE;
      }
    }

    return t;
  }

  private static void writeTag(DataOutputStream writer, Object tag, int tagsStartBytes,
      Map<String, Integer> tagOffsetBytes, int intTagsStartBytes, Map<Integer, Integer> intTagOffsetBytes,
      int doubleTagsStartBytes, Map<Double, Integer> doubleTagOffsetBytes) throws IOException {
    TagType t = getTagType(tag);

    writer.write(TagType.byteRep(t));

    switch (t) {
    case INT:
      if (tag instanceof Integer) {
        writer.writeInt(intTagsStartBytes + intTagOffsetBytes.get((int) tag));
      } else {
        writer.writeInt(intTagsStartBytes + intTagOffsetBytes.get(Integer.parseInt(tag.toString())));
      }
      break;
    case DOUBLE:
      if (tag instanceof Double) {
        writer.writeInt(intTagsStartBytes + doubleTagOffsetBytes.get((double) tag));
      } else {
        writer.writeInt(doubleTagsStartBytes + doubleTagOffsetBytes.get(Double.parseDouble(tag.toString())));
      }
      break;
    default:
      writer.writeInt(tagsStartBytes + tagOffsetBytes.get(tag.toString()));
      break;
    }
  }

  private static int elementSizeBytes(GenomicElement g) {
    // Size of entry in bytes
    // Prefix each read with 'b' for ease of finding
    int s = 0; // Geb.INT_BYTES;

    // id
    // s += GEBReader.INT_BYTES;

    // type reference
    s += GEBReader.INT_BYTES;

    // chr (String) start end
    s += GEBReader.INT_BYTES + GEBReader.INT_BYTES + GEBReader.INT_BYTES;

    // Strand
    s += 1;

    // number of exons + exon starts and ends
    // s += 1 + e.getChildCount(GenomicEntity.EXON) * Geb.INT_BYTES * 2;

    // 1 byte for number of ids and each id needs 1 byte for type and 2 ints
    // (key, value) location
    s += 1 + Math.min(g.getPropertyCount(), GEBReader.MAX_TAGS) * (1 + 2 * GEBReader.INT_BYTES);

    /*
     * for (Entry<String, String> id : g.getIds()) { // 1 byte for length of type
     * then that number of bytes space + // 1 byte for length of value + that number
     * of bytes
     * 
     * // Address of string s += Geb.INT_BYTES; //varcharSize(id.getKey()) +
     * varcharSize(id.getValue()); }
     */

    // number of tags (must be fewer than 256)
    s += 1 + Math.min(g.getTagCount(), GEBReader.MAX_TAGS) * (1 + GEBReader.INT_BYTES);

    /*
     * for (String tag : g.getTags()) { // 1 byte for length of tag plus that number
     * of bytes to store // tag value
     * 
     * // Address of string s += Geb.INT_BYTES; //varcharSize(tag); }
     */

    // Number of children (short)
    s += 2;

    int size = 0;

    for (GenomicType childType : g.getChildTypes()) {
      for (GenomicElement child : g.getChildren(childType)) {
        s += elementSizeBytes(child);

        if (size++ == GEBReader.MAX_CHILDREN) {
          break;
        }
      }
    }

    return s;
  }

  private static <T extends GenomicElement> void tagSizesBytes(final Iterable<T> features,
      IterMap<String, Integer> sizeMap, IterMap<Integer, Integer> intSizeMap, IterMap<Double, Integer> doubleSizeMap,
      IterMap<String, Integer> offsetMap, IterMap<Integer, Integer> intOffsetMap,
      IterMap<Double, Integer> doubleOffsetMap) {

    for (GenomicElement feature : features) {
      tagSizesBytes(feature, sizeMap, intSizeMap, doubleSizeMap);
    }

    int n = 0;

    for (Entry<String, Integer> item : sizeMap) {
      String s = item.getKey();

      if (!TextUtils.isInt(s) && !TextUtils.isDouble(s)) {
        offsetMap.put(item.getKey(), n);
        n += item.getValue();
      }
    }

    n = 0;

    for (Entry<Integer, Integer> e : intSizeMap) {
      intOffsetMap.put(e.getKey(), n);
      n += e.getValue();
    }

    n = 0;

    for (Entry<Double, Integer> e : doubleSizeMap) {
      doubleOffsetMap.put(e.getKey(), n);
      n += e.getValue();
    }
  }

  private static void tagSizesBytes(GenomicElement feature, Map<String, Integer> sizeMap,
      Map<Integer, Integer> intSizeMap, Map<Double, Integer> doubleSizeMap) {

    Deque<GenomicElement> stack = new ArrayDeque<GenomicElement>();

    stack.push(feature);

    while (!stack.isEmpty()) {
      GenomicElement f = stack.pop();

      String v = f.getType().toString();
      sizeMap.put(v, varcharSize(v));

      v = f.getChr().toString();
      sizeMap.put(v, varcharSize(v));

      for (String k : f.getPropertyNames()) {
        sizeMap.put(k, varcharSize(k));

        tagSizeBytes(f.getProperty(k), sizeMap, intSizeMap, doubleSizeMap);
      }

      for (String t : f.getTags()) {
        tagSizeBytes(t, sizeMap, intSizeMap, doubleSizeMap);
      }

      for (GenomicType t : f.getChildTypes()) {
        for (GenomicElement c : f.getChildren(t)) {
          // tagSizesBytes(c, sizeMap, intSizeMap, doubleSizeMap);
          stack.push(c);
        }
      }
    }
  }

  private static void tagSizeBytes(Object tag, Map<String, Integer> sizeMap, Map<Integer, Integer> intSizeMap,
      Map<Double, Integer> doubleSizeMap) {

    if (tag instanceof Integer) {
      intSizeMap.put((int) tag, GEBReader.INT_BYTES);
    } else if (tag instanceof Double) {
      doubleSizeMap.put((double) tag, GEBReader.DOUBLE_BYTES);
    } else {
      String v = tag.toString();

      if (TextUtils.isInt(v)) {
        intSizeMap.put(Integer.parseInt(v), GEBReader.INT_BYTES);
      } else if (TextUtils.isDouble(v)) {
        doubleSizeMap.put(Double.parseDouble(v), GEBReader.DOUBLE_BYTES);
      } else {
        sizeMap.put(v, varcharSize(v));
      }
    }
  }

  private void writeData(IterTreeMap<String, Integer> tagOffsetBytesMap,
      IterTreeMap<Integer, Integer> intTagOffsetBytesMap, IterTreeMap<Double, Integer> doubleOffsetBytesMap)
      throws IOException {
    Path file = mDir.resolve(DataReader.getFileName(mPrefix)); // PathUtils.getPath(mGenome

    DataOutputStream writer = FileUtils.newDataOutputStream(file);

    LOG.info("Writing data to {}...", file);

    writer.writeInt(GEBReader.CHECK);
    writer.writeByte(GEBReader.VERSION);
    writer.writeInt(mWindow);

    // Write Strings first
    writer.writeInt(tagOffsetBytesMap.size());

    for (Entry<String, Integer> e : tagOffsetBytesMap) {
      // Each gene is given an id to make it easier to check which genes
      // we have found when decoding

      writeVarchar(e.getKey(), writer);
    }

    // Write ints
    writer.writeInt(intTagOffsetBytesMap.size());

    for (Entry<Integer, Integer> e : intTagOffsetBytesMap) {
      // Each gene is given an id to make it easier to check which genes
      // we have found when decoding
      writer.writeInt(e.getKey());
    }

    // Write double
    writer.writeInt(doubleOffsetBytesMap.size());

    for (Entry<Double, Integer> e : doubleOffsetBytesMap) {
      // Each gene is given an id to make it easier to check which genes
      // we have found when decoding
      writer.writeDouble(e.getKey());
    }

    writer.close();
  }

  private void writeIndex() throws IOException {
    Path file = mDir.resolve(GEBReader.getIndexFileName(mPrefix)); // PathUtils.getPath(mGenome

    JsonObject root = new JsonObject();

    root.add("name", mPrefix);
    root.add("window", mWindow);

    JsonObject go = new JsonObject();
    go.add("name", "Human");
    go.add("build", "hg19");
    root.add("genome", go);

    Json.prettyWrite(root, file);

    // BufferedWriter writer = FileUtils.newBufferedWriter(file);

    // writer.write(mPrefix);
    // writer.newLine();

    // writer.close();
  }

  private void writeBins(Chromosome chr, IterMap<Integer, List<GenomicElement>> binsMap,
      IterMap<Integer, Integer> binSizeBytes, IterMap<GenomicElement, Integer> elementOffsetBytes) throws IOException {
    Path file = mDir.resolve(BinReader.getFileName(mPrefix, chr)); // PathUtils.getPath(mGenome

    DataOutputStream writer = FileUtils.newDataOutputStream(file);

    LOG.info("Writing bins to {}...", file);

    int n = binsMap.size();
    int minBin = binsMap.keySet().iterator().next();

    // Version
    // The first int should be 42 so that you can tell whether the
    // endian is correct
    writer.writeInt(GEBReader.CHECK);
    writer.writeByte(GEBReader.VERSION);
    writer.writeInt(mWindow);
    writer.writeInt(minBin);
    writer.writeInt(n);

    int binAddressesWidthBytes = n * GEBReader.INT_BYTES;

    int offset = BinReader.HEADER_BYTES_OFFSET + binAddressesWidthBytes;

    for (Entry<Integer, Integer> item : binSizeBytes) {
      // Write address to bin containing element ids
      writer.writeInt(offset);

      offset += item.getValue();
    }

    // Write out each bin
    for (Entry<Integer, Integer> item : binSizeBytes) {
      List<GenomicElement> bins = binsMap.get(item.getKey());
      // writer.write(GEBReader.BLOCK_SEPARATOR);

      // Number of gene addresses in the bin
      writer.writeInt(bins.size());

      // Write the addresses to each gene
      for (GenomicElement e : bins) {
        writer.writeInt(ElementReader.HEADER_BYTES_OFFSET + elementOffsetBytes.get(e));
      }
    }

    writer.close();
  }

  private void writeBTree(Chromosome chr, IterMap<Integer, List<GenomicElement>> binsMap,
      IterMap<Integer, Integer> binSizeBytes, IterMap<GenomicElement, Integer> elementOffsetBytes) throws IOException {

    // bins map contains just the in use bins

    // list of the bins actually in use
    List<Integer> bins = CollectionUtils.sort(binsMap.keySet());

    int me = bins.size() - 1;

    IterMap<Integer, Integer> binOffsetBytes = new IterTreeMap<Integer, Integer>();

    int offset = 0;

    for (Entry<Integer, List<GenomicElement>> item : binsMap) {
      binOffsetBytes.put(item.getKey(), offset);

      offset += GEBReader.INT_BYTES; // binSizeBytes.get(item.getKey());
    }

    Deque<BTreeNode<GenomicElement>> q = new ArrayDeque<BTreeNode<GenomicElement>>();
    Deque<Integer> indexStack = new ArrayDeque<Integer>();

    BTreeNode<GenomicElement> root = new BTreeNode<GenomicElement>(me / 2);
    q.push(root);
    indexStack.push(me);
    indexStack.push(0);

    int treeSizeOffset = 0; // BTreeReader.HEADER_BYTES_OFFSET;

    IterMap<BTreeNode<GenomicElement>, Integer> nodeOffsetBytes = new IterHashMap<BTreeNode<GenomicElement>, Integer>();

    while (!q.isEmpty()) {
      BTreeNode<GenomicElement> node = q.pop();
      int s = indexStack.pop();
      int e = indexStack.pop();

      int w = btreeNodeSizeBytes(node);

      nodeOffsetBytes.put(node, treeSizeOffset);

      treeSizeOffset += w;

      // Add objects to node
      for (GenomicElement be : binsMap.get(bins.get(node.getIndex()))) {
        node.add(be);
      }

      // System.err.println("node " + node.getIndex());

      if (node.getIndex() < e) {
        int ni = node.getIndex() + 1;
        BTreeNode<GenomicElement> child = new BTreeNode<GenomicElement>((ni + e) / 2);

        node.setC2(child);

        q.push(child);
        indexStack.push(e);
        indexStack.push(ni);

        // System.err.println("> " + node.getIndex() + " " + child.getIndex() +
        // " "
        // + ni + " " + e);
      }

      // if we the same as s, then there is nothing to be done
      if (node.getIndex() > s) {
        int ni = node.getIndex() - 1;

        BTreeNode<GenomicElement> child = new BTreeNode<GenomicElement>((s + ni) / 2);

        node.setC1(child);

        q.push(child);
        indexStack.push(ni);
        indexStack.push(s);

        // System.err.println("< " + node.getIndex() + " " + child.getIndex() +
        // " "
        // + s + " " + ni);
      }
    }

    //
    // Write to file
    //

    Path file = mDir.resolve(BTreeReader.getFileName(mPrefix, chr));

    LOG.info("Writing btree to {}...", file);

    DataOutputStream writer = FileUtils.newDataOutputStream(file);

    writer.writeInt(GEBReader.CHECK);
    writer.writeByte(GEBReader.VERSION);
    writer.writeInt(mWindow);
    writer.writeInt(BTreeReader.HEADER_BYTES_OFFSET + treeSizeOffset);

    q = new ArrayDeque<BTreeNode<GenomicElement>>();

    q.push(root);

    while (!q.isEmpty()) {
      BTreeNode<GenomicElement> node = q.pop();

      int bin = bins.get(node.getIndex());

      writer.writeInt(bin);

      // Write address to bin
      writer.writeInt(BTreeReader.HEADER_BYTES_OFFSET + treeSizeOffset + binOffsetBytes.get(bin));

      if (node.getC1() != null) {
        writer.writeInt(BTreeReader.HEADER_BYTES_OFFSET + nodeOffsetBytes.get(node.getC1()));
      } else {
        writer.writeInt(0);
      }

      if (node.getC2() != null) {
        writer.writeInt(BTreeReader.HEADER_BYTES_OFFSET + nodeOffsetBytes.get(node.getC2()));
      } else {
        writer.writeInt(0);
      }

      // Process the children
      if (node.getC2() != null) {
        q.push(node.getC2());
      }

      if (node.getC1() != null) {
        q.push(node.getC1());
      }
    }

    int binAddressesOffset = binsMap.size() * GEBReader.INT_BYTES;

    offset = BTreeReader.HEADER_BYTES_OFFSET + treeSizeOffset + binAddressesOffset;

    for (Entry<Integer, List<GenomicElement>> item : binsMap) {
      writer.writeInt(offset);

      offset += binSizeBytes.get(item.getKey());
    }

    for (Entry<Integer, List<GenomicElement>> item : binsMap) {
      // Number of gene addresses in the bin
      writer.writeInt(item.getValue().size());

      // Write the addresses to each gene
      for (GenomicElement e : item.getValue()) {
        writer.writeInt(ElementReader.HEADER_BYTES_OFFSET + elementOffsetBytes.get(e));
      }
    }

    writer.close();
  }

  private static final int btreeNodeSizeBytes(BTreeNode<GenomicElement> node) {
    // Bin
    int s = GEBReader.INT_BYTES;

    // offset in bin array
    s += GEBReader.INT_BYTES;

    // Addresses of two children
    s += BTreeReader.BTREE_CHILD_ADDRESSES_BYTES; // node.getChildCount() *
    // BTreeReader.BTREE_TREE_PREFIX_BYTES;

    // address of bin
    // s += GEBReader.INT_BYTES * (1 + node.getObjectCount());

    return s;
  }

  private <T extends GenomicElement> void writeRadix(IterMap<Chromosome, Set<T>> elements,
      IterMap<GenomicElement, Integer> elementOffsetBytes) throws IOException {

    // The header + the space occupied by the bin addresses + the space
    // occupied by the bins
    // int treeOffset = GFBGenes.RADIX_BYTES_OFFSET;

    RadixNode<GenomicElement> root = new RadixNode<GenomicElement>();

    Deque<RadixNode<GenomicElement>> q;

    for (Entry<Chromosome, Set<T>> c : elements) {
      Set<T> features = c.getValue(); // chrMap.get(chr);

      // Make a radix tree

      for (GenomicElement e : features) {
        root.add(e.getChr().toString(), e);
        // root.add(Integer.toString(e.getStart()), e);
        // root.add(Integer.toString(e.getEnd()), e);

        // Index properties
        for (String item : e.getPropertyNames()) {
          String v = e.getProperty(item);

          root.add(v, e);
        }

        // Tags
        for (String tag : e.getTags()) {
          root.add(tag.toString(), e);
        }
      }
    }

    // Calculate space occupied by tree

    IterMap<RadixNode<GenomicElement>, Integer> nodeOffsetBytes = new IterHashMap<RadixNode<GenomicElement>, Integer>();

    // Record offsets in the order we encounter nodes since mutliple nodes
    // may have same
    // List<Integer> nodeOffsetBytes = new ArrayList<Integer>();

    q = new ArrayDeque<RadixNode<GenomicElement>>();

    q.push(root);

    int offset = 0;

    while (!q.isEmpty()) {
      RadixNode<GenomicElement> node = q.pop();

      // Each node consists of a char + the number of children +
      // address to each child
      int w = radixNodeSizeBytes(node);

      nodeOffsetBytes.put(node, offset);

      offset += w;

      // Push all the children on
      for (Entry<Character, RadixNode<GenomicElement>> item : node.getChildren()) {
        q.push(item.getValue());
      }
    }

    //
    // Write to file
    //

    Path file = mDir.resolve(RadixReader.getFileName(mPrefix));

    LOG.info("Writing radix to {}...", file);

    DataOutputStream writer = FileUtils.newDataOutputStream(file);

    writer.writeInt(GEBReader.CHECK);
    writer.writeByte(GEBReader.VERSION);
    writer.writeInt(mWindow);

    q = new ArrayDeque<RadixNode<GenomicElement>>();

    q.push(root);

    while (!q.isEmpty()) {
      RadixNode<GenomicElement> node = q.pop();

      // The number of children
      // writer.writeInt(node.getChildCount());
      writer.writeByte(node.getChildCount());

      // System.err.println("write " + node.getChar() + " " +
      // node.getChildNames());

      for (Entry<Character, RadixNode<GenomicElement>> item : node.getChildren()) {
        char c = item.getKey();

        // The char this represents
        writer.write(c);

        writer.writeInt(RadixReader.HEADER_BYTES_OFFSET + nodeOffsetBytes.get(item.getValue()));
      }

      // Write addresses of elements

      writer.writeInt(node.getExactObjects().size());

      for (GenomicElement e : node.getExactObjects()) {
        writer.writeInt(ElementReader.HEADER_BYTES_OFFSET + elementOffsetBytes.get(e));
      }

      writer.writeInt(node.getObjects().size());

      for (GenomicElement e : node.getObjects()) {
        writer.writeInt(ElementReader.HEADER_BYTES_OFFSET + elementOffsetBytes.get(e));
      }

      // Process the children
      for (Entry<Character, RadixNode<GenomicElement>> item : node.getChildren()) {
        q.push(item.getValue());
      }
    }

    writer.close();
  }

  private static final int radixNodeSizeBytes(RadixNode<GenomicElement> node) {
    // Number of children
    int s = 1;

    // n x (1 byte char + int address)
    s += node.getChildCount() * RadixReader.RADIX_TREE_PREFIX_BYTES;

    // exact matches
    s += GEBReader.INT_BYTES * (1 + node.getExactObjects().size());

    // partial matches
    s += GEBReader.INT_BYTES * (1 + node.getObjects().size());

    return s;
  }
}
