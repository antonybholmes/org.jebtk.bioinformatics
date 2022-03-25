/**
 * Copyright 2018 Antony Holmes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.jebtk.bioinformatics.genomic;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jebtk.core.collections.ArrayUtils;
import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.Join;
import org.jebtk.core.text.TextUtils;

/**
 * Encode genes in binary binned format.
 *
 * @author Antony Holmes
 */
public class GFBGenes extends SingleGenesDB {

  private static final long serialVersionUID = 1L;

  /** The buffer. */
  private final byte[] mBuffer = new byte[256];

  /** The bin addresses. */
  private final int[] mBinAddresses = new int[300000];

  /** The gene addresses. */
  private final int[] mGeneAddresses = new int[200000];

  /** keep track of which addresses are used during a search **/
  // private final boolean[] mUsed = new boolean[200000];

  private final GenomicElement[] mGenes = new GenomicElement[200000];

  /**
   * The Constant INT_BYTES represents the bytes used by a 32bit number (8 * 4)
   */
  public static final int INT_BYTES = 4;

  /** The Constant VERSION_OFFSET. */
  private static final int VERSION_BYTE_OFFSET = INT_BYTES;

  /** The Constant DESCRIPTION_OFFSET. */
  // private static final int DESCRIPTION_BYTE_OFFSET = VERSION_BYTE_OFFSET + 1;

  /** The Constant DESCRIPTION_BYTES. */
  // public static final int DESCRIPTION_BYTES = 16;

  /** The Constant DESCRIPTION_BYTES_USABLE. */
  // public static final int DESCRIPTION_BYTES_USABLE = DESCRIPTION_BYTES - 1;

  /** The Constant WINDOW_BYTE_OFFSET. */
  public static final int GENES_BYTES_OFFSET = VERSION_BYTE_OFFSET + 1;

  public static final int RADIX_BYTES_OFFSET = GENES_BYTES_OFFSET + INT_BYTES;

  /** The Constant WINDOW_BYTE_OFFSET. */
  private static final int WINDOW_BYTE_OFFSET = RADIX_BYTES_OFFSET;

  /** The Constant BINS_BYTE_OFFSET. */
  private static final int BINS_BYTE_OFFSET = WINDOW_BYTE_OFFSET + INT_BYTES;

  /** The Constant HEADER_BYTES is where the bins start */
  public static final int HEADER_BYTES_OFFSET = BINS_BYTE_OFFSET + INT_BYTES; // GFBGenes.INT_BYTES
  // + 1 +
  // GFBGenes.DESCRIPTION_BYTES
  // +
  // GFBGenes.INT_BYTES
  // +
  // GFBGenes.INT_BYTES;

  public static final int RADIX_TREE_PREFIX_BYTES = 1 + INT_BYTES;

  public static final int GENOMIC_TYPE_GENE = 1;
  public static final int GENOMIC_TYPE_TRANSCRIPT = 2;
  public static final int GENOMIC_TYPE_EXON = 4;

  private final Path mDir;

  /** The m window. */
  private int mWindow;

  private int mAddress;

  private long mPos;

  private int mAddress2;

  /**
   * Instantiates a new GFB genes.
   *
   * @param genome the genome
   * @param window the window
   * @param dir    the dir
   */
  public GFBGenes(Genome genome, int window, Path dir) {
    super(genome);

    mWindow = window;
    mDir = dir;
  }

  /**
   * Find genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public List<GenomicElement> findGenes(String region, int minBp) throws IOException {
    return findGenes(region, GenomicType.TRANSCRIPT, minBp);
  }

  public List<GenomicElement> findGenes(String region, GenomicType type, int minBp) throws IOException {
    return find(null, GenomicRegion.parse(mGenome, region), type, minBp);
  }

  /**
   * Find genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  @Override
  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
    return findGenes(region.mChr, region.mStart, region.mEnd, type);
  }

  /**
   * Find genes.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public List<GenomicElement> findGenes(Chromosome chr, int start, int end, GenomicType type) {

    int n = _findGenes(chr, start, end, type);

    return overlappingGenes(chr, start, end, n);
  }

  /**
   * Find genes in the blocks spanning the coordinates. These are the genes most
   * likely to be overlapping the region of interest. A further test is required
   * to test for overlap. This method is designed to narrow down the list of genes
   * 
   * @param chr
   * @param start
   * @param end
   * @param genes
   * @return
   * @throws IOException
   */
  private int _findGenes(Chromosome chr, int start, int end, GenomicType type) {
    int sb = start / mWindow;
    int eb = end / mWindow;

    int n = 0;

    int[] bins = ArrayUtils.array(sb, eb);

    Path file = mDir.resolve(getFileName(mGenome, chr, mWindow));

    RandomAccessFile reader;

    try {
      reader = FileUtils.newRandomAccess(file);

      try {
        n = binAddressesFromBins(reader, bins);

        n = geneAddressesFromBins(reader, n);

        n = genesFromGeneAddresses(reader, n, type);

      } finally {
        reader.close();
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    return n;
  }

  @Override
  public List<GenomicElement> getElements() {
    Path file = mDir.resolve(getRadixFileName(mGenome));

    int n = 0;

    RandomAccessFile reader;
    try {
      reader = FileUtils.newRandomAccess(file);

      try {
        int address = readGenesAddress(reader);

        n = readAllGenes(reader, address, GenomicType.GENE);
      } finally {
        reader.close();
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    return ArrayUtils.toList(mGenes, n);
  }

  public List<GenomicElement> getGenes(String id) throws IOException {
    return getGenes(id, GenomicType.TRANSCRIPT);
  }

  @Override
  public List<GenomicElement> getElements(Genome genome, String search, GenomicType type) {
    return getGenes(search, GenomicType.TRANSCRIPT);
  }

  public List<GenomicElement> getGenes(String search, GenomicType type) {
    Path file = mDir.resolve(getRadixFileName(mGenome));

    int n = 0;

    RandomAccessFile reader;
    try {
      reader = FileUtils.newRandomAccess(file);

      try {
        n = geneAddressesFromRadix(reader, search);

        n = genesFromGeneAddresses(reader, n, type);
      } finally {
        reader.close();
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    return ArrayUtils.toList(mGenes, n);
  }

  @Override
  public List<GenomicElement> closest(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    return findClosestGenes(region, type, minBp);
  }

  public List<GenomicElement> findClosestGenes(GenomicRegion region, GenomicType type, int minBp) throws IOException {
    return findClosestGenes(region.mChr, region.mStart, region.mEnd, type, minBp);
  }

  public List<GenomicElement> findClosestGenes(Chromosome chr, int start, int end, int minBp) throws IOException {
    return findClosestGenes(chr, start, end, GenomicType.TRANSCRIPT, minBp);
  }

  public List<GenomicElement> findClosestGenes(Chromosome chr, int start, int end, GenomicType type, int minBp)
      throws IOException {

    int n = _findGenes(chr, start, end, type);

    return findClosestGenes(chr, start, end, n, minBp);
  }

  private List<GenomicElement> findClosestGenes(Chromosome chr, int start, int end, int n, int minBp) {
    if (n < 1) {
      return Collections.emptyList();
    }

    List<GenomicElement> ret = new ArrayList<GenomicElement>(n);

    int mid = GenomicRegion.mid(start, end);
    int minD = Integer.MAX_VALUE;

    for (int i = 0; i < n; ++i) {
      GenomicElement gene = mGenes[i];

      int d = Math.abs(mid - gene.getStart()); // GenomicRegion.mid(gene));

      if (d < minD) {
        minD = d;
      }
    }

    for (int i = 0; i < n; ++i) {
      GenomicElement gene = mGenes[i];

      int d = Math.abs(mid - gene.getStart()); // GenomicRegion.mid(gene));

      if (d == minD) {
        ret.add(gene);
      }
    }

    return ret;
  }

  /**
   * Gets the overlapping genes.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @param genes the genes
   * @return the overlapping genes
   */
  private List<GenomicElement> overlappingGenes(Chromosome chr, int start, int end, int n) {
    return overlappingGenes(chr, start, end, n, 1);
  }

  /**
   * Gets the overlapping genes.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @param genes the genes
   * @param minBp the min bp
   * @return the overlapping genes
   */
  private List<GenomicElement> overlappingGenes(Chromosome chr, int start, int end, int n, int minBp) {
    if (n < 1) {
      return Collections.emptyList();
    }

    List<GenomicElement> ret = new ArrayList<GenomicElement>(n);

    for (int i = 0; i < n; ++i) {
      GenomicElement gene = mGenes[i];

      GenomicRegion overlap = GenomicRegion.overlap(chr, start, end, gene);

      if (overlap != null && (overlap.getLength() >= minBp)) {
        ret.add(gene);
      }
    }

    return ret;
  }

  /**
   * Gets the bin addresses in the files for a list of bins.
   *
   * @param reader    the reader
   * @param bins      a list of bins.
   * @param addresses array that will be populated with the bin addresses in the
   *                  same order as bins.
   * @return the number of bins.
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int binAddressesFromBins(RandomAccessFile reader, final int[] bins) throws IOException {
    int n = bins.length;

    for (int i = 0; i < n; ++i) {
      mBinAddresses[i] = readBinAddress(reader, bins[i]);
    }

    return n;
  }

  /**
   * Get the gene addresses from a selection of bin addresses, removing
   * duplicates. A Gene address corresponds to a gene transcript.
   *
   * @param reader        a GFB binary file.
   * @param binAddresses  an array of bin addresses
   * @param n             how many bin addresses are in use.
   * @param geneAddresses an array of gene addresses that will be populated with
   *                      results.
   * @return the gene addresses
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int geneAddressesFromBins(RandomAccessFile reader, int n) throws IOException {

    int retn = 0;

    // bin address
    int ba;

    // Gene address
    int ga;

    // How many genes are in the bin
    int size;

    Set<Integer> used = new HashSet<Integer>();

    for (int i = 0; i < n; ++i) {
      ba = mBinAddresses[i];

      reader.seek(ba);

      size = reader.readInt();

      // Read how many addresses are in the bin and then extract them
      for (int j = 0; j < size; ++j) {
        ga = reader.readInt();

        // Remove duplicates
        if (!used.contains(ga)) {
          mGeneAddresses[retn++] = ga;
          used.add(ga);
        }
      }
    }

    // Return the number of genes we found
    return retn;
  }

  /**
   * Loads genes from an array of gene addresses
   *
   * @param reader    the reader
   * @param chr       the chr
   * @param addresses the addresses
   * @param n         the n
   * @return the genes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int genesFromGeneAddresses(RandomAccessFile reader, int n, GenomicType type) throws IOException {
    int retn = 0;

    // Clear the used array
    // Arrays.fill(mUsed, false);

    for (int i = 0; i < n; ++i) {
      // Skip to address of gene
      reader.seek(mGeneAddresses[i]);

      retn = readGene(reader, retn, type);
    }

    return retn;
  }

  /**
   * Read all the genes.
   * 
   * @param reader
   * @param address
   * @param genes
   * @return
   * @throws IOException
   */
  private int readAllGenes(RandomAccessFile reader, final int address, GenomicType type) throws IOException {
    int retn = 0;

    // Clear the used array
    // Arrays.fill(mUsed, false);

    reader.seek(address);

    int n = reader.readInt();

    for (int i = 0; i < n; ++i) {
      // Read all the genes
      retn = readGene(reader, retn, type);
    }

    return retn;
  }

  private int readGene(RandomAccessFile reader, int index, GenomicType type) throws IOException {
    // System.err.println("Reading gene");

    GenomicElement gene = readEntity(reader, GenomicType.GENE);

    int n = reader.read();

    for (int i = 0; i < n; ++i) {
      // If type requested is not a gene, pass null to indicate that the
      // transcripts should not add themselves to the gene
      index = readTranscript(reader, index, type, type == GenomicType.GENE ? gene : null);
    }

    if (type == GenomicType.GENE) {
      mGenes[index++] = gene;
    }

    return index;
  }

  private int readTranscript(RandomAccessFile reader, int index, GenomicType type, GenomicElement gene)
      throws IOException {
    // System.err.println("Reading transcript");

    GenomicElement transcript = readEntity(reader, GenomicType.TRANSCRIPT);

    int n = reader.read();

    for (int i = 0; i < n; ++i) {
      index = readExon(reader, index, type, gene, type != GenomicType.EXON ? transcript : null);
    }

    if (gene != null) {
      gene.addChild(transcript);
    }

    if (type == GenomicType.TRANSCRIPT) {
      mGenes[index++] = transcript;
    }

    return index;
  }

  private int readExon(RandomAccessFile reader, int index, GenomicType type, GenomicElement gene,
      GenomicElement transcript) throws IOException {
    // System.err.println("Reading exon");

    GenomicElement exon = readEntity(reader, GenomicType.EXON);

    // skip byte for exon child count, since exons cannot have children
    // so this byte is always set to zero.
    reader.skipBytes(1);

    if (transcript != null) {
      transcript.addChild(exon);
    }

    if (type == GenomicType.EXON) {
      mGenes[index++] = exon;
    }

    return index;
  }

  private GenomicElement readEntity(RandomAccessFile reader, GenomicType type) throws IOException {

    // Skip id (int) and type (byte)
    reader.skipBytes(INT_BYTES + 1); // .readInt();

    GenomicRegion l = readLocation(reader);

    Strand strand = readStrand(reader);

    GenomicElement gene;

    switch (type) {
    case GENE:
      gene = new Gene(l, strand);
      break;
    case EXON:
      gene = new Exon(l, strand);
      break;
    default:
      // Assume everything is a transcript
      gene = new Transcript(l, strand);
      break;
    }

    readIds(reader, gene);

    readTags(reader, gene);

    return gene;
  }

  private GenomicRegion readLocation(RandomAccessFile reader) throws IOException {
    Chromosome chr = Chromosome.newChr(readVarchar(reader));

    int start = reader.readInt();
    int end = reader.readInt();

    return GenomicRegion.create(chr, start, end);
  }

  private Strand readStrand(RandomAccessFile reader) throws IOException {
    int strand = reader.read();

    return getStrand(strand);
  }

  /**
   * Read ids and add them to a genomic entity, e.g. read gene symbol and add to a
   * transcript.
   * 
   * @param reader
   * @param e
   * @return
   * @throws IOException
   */
  private int readIds(RandomAccessFile reader, GenomicElement e) throws IOException {
    int n = reader.readByte();

    for (int i = 0; i < n; ++i) {
      readId(reader, e);
    }

    return n;
  }

  private void readId(RandomAccessFile reader, GenomicElement e) throws IOException {

    mAddress = reader.readInt();
    mAddress2 = reader.readInt();

    mPos = reader.getFilePointer();

    e.setProperty(readString(reader, mAddress), readString(reader, mAddress2));

    reader.seek(mPos);
  }

  /**
   * Load tags associated with an entity.
   * 
   * @param reader
   * @param e
   * @return
   * @throws IOException
   */
  private int readTags(RandomAccessFile reader, GenomicElement e) throws IOException {
    int n = reader.readByte();

    // System.err.println("Read tags " + n);

    for (int i = 0; i < n; ++i) {
      readTag(reader, e);
    }

    return n;
  }

  private void readTag(RandomAccessFile reader, GenomicElement e) throws IOException {
    e.addTag(readTag(reader));
  }

  private String readTag(RandomAccessFile reader) throws IOException {

    mAddress = reader.readInt();

    mPos = reader.getFilePointer();

    String ret = readString(reader, mAddress);

    reader.seek(mPos);

    return ret;
  }

  private String readString(RandomAccessFile reader, int address) throws IOException {

    reader.seek(address);

    return readVarchar(reader);
  }

  /**
   * Read a variable number of bytes to create a GenomicEntity.
   *
   * @param reader the reader
   * @return the string
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private String readVarchar(RandomAccessFile reader) throws IOException {
    // System.err.println("vchar " + reader.getFilePointer());

    // First int tells us the length of the string
    int n = reader.readByte();

    // Read n bytes into the buffer
    reader.read(mBuffer, 0, n);

    // Create string from buffer
    String s = new String(mBuffer, 0, n, StandardCharsets.UTF_8);

    return s;
  }

  /*
   * private static int getRadixRootGeneAddresses(RandomAccessFile reader, int[]
   * geneAddresses) throws IOException { reader.seek(RADIX_BYTES_OFFSET);
   * 
   * // Number of children int n = reader.readInt();
   * 
   * // skip tree addresses // reader.seek(RADIX_HEADER_BYTES + N_BYTES + n *
   * RADIX_TREE_PREFIX_BYTES); reader.skipBytes(n * RADIX_TREE_PREFIX_BYTES);
   * 
   * n = reader.readInt();
   * 
   * for (int i = 0; i < n; ++i) { geneAddresses[i] = reader.readInt(); }
   * 
   * // Return the number of genes we found return n; }
   */

  /**
   * Populates internal address list with the gene addresses from the radix tree
   * so that a text search can be performed.
   * 
   * @param reader
   * @param id
   * @return
   * @throws IOException
   */
  private int geneAddressesFromRadix(RandomAccessFile reader, String id) throws IOException {

    char[] chars = id.toLowerCase().toCharArray();

    // Find the tree start
    reader.seek(RADIX_BYTES_OFFSET);

    char leafc;
    int address = 0;
    int n;

    boolean found = false;

    for (char c : chars) {

      // Number of children
      n = reader.read(); // .readInt();

      // System.err.println("c " + c + " " + n);

      // assume we won't find a match
      found = false;

      for (int i = 0; i < n; ++i) {
        leafc = (char) reader.read(); // readByte();
        address = reader.readInt();

        if (leafc == c) {
          // we did find a match so keep going
          found = true;
          reader.seek(address);
          break;
        }
      }

      if (!found) {
        break;
      }
    }

    if (!found) {
      return 0;
    }

    // This means we kept finding a child matching the prefix so the seek
    // is at the beginning of a node either because we ran out of chars
    // or nodes. In this case we must skip over the child addresses and
    // just look at the addresses of the genes

    // skip past number of addresses and the addresses themselves
    // reader.seek(address + N_BYTES + reader.readInt() *
    // RADIX_TREE_PREFIX_BYTES);
    reader.skipBytes(reader.read() * RADIX_TREE_PREFIX_BYTES); // readInt()

    // Should be on a node that is hopefully matches our search term
    // Since we checked all the children, the seek is at the position of
    // the gene addressses so we can get them

    n = reader.readInt();

    for (int i = 0; i < n; ++i) {
      mGeneAddresses[i] = reader.readInt();
    }

    // Return the number of genes we found
    return n;
  }

  @Override
  public Iterable<String> getNames() throws IOException {
    List<GenomicElement> genes = getElements();

    List<String> ret = new UniqueArrayList<String>(genes.size());

    for (GenomicElement gene : genes) {
      ret.add(((GenomicEntity) gene).getSymbol());
    }

    return ret;
  }

  @Override
  public void add(GenomicElement element) {
    // Do nothing
  }

  /**
   * Gets the file name.
   *
   * @param genome the genome
   * @param chr    the chr
   * @param window the window
   * @return the file name
   */
  public static final Path getFileName(Genome genome, Chromosome chr, int window) {
    return PathUtils
        .getPath(Join.on('.').values(genome.getAssembly(), chr, TextUtils.cat("w", window), "gfb").toString());
  }

  public static final Path getRadixFileName(Genome genome) {
    return PathUtils.getPath(genome.getAssembly() + ".rgfb");
  }

  public static int readCheckNum(RandomAccessFile reader) throws IOException {
    reader.seek(0);

    return reader.readByte();
  }

  /**
   * Gets the version.
   *
   * @param reader the reader
   * @return the version
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static int readVersion(RandomAccessFile reader) throws IOException {
    reader.seek(VERSION_BYTE_OFFSET);

    return reader.readByte();
  }

  /*
   * public String getDescription(RandomAccessFile reader) throws IOException {
   * reader.seek(DESCRIPTION_BYTE_OFFSET); int n = reader.readByte();
   * reader.read(mBuffer, 0, DESCRIPTION_BYTES_USABLE); return new String(mBuffer,
   * 0, n); }
   */

  /**
   * Gets the bins.
   *
   * @param reader the reader
   * @return the bins
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static int readBinCount(RandomAccessFile reader) throws IOException {
    reader.seek(BINS_BYTE_OFFSET);

    return reader.readInt();
  }

  /**
   * Gets the window.
   *
   * @param reader the reader
   * @return the window
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static int readWindow(RandomAccessFile reader) throws IOException {
    reader.seek(WINDOW_BYTE_OFFSET);

    return reader.readInt();
  }

  /**
   * Gets the address of where the genes begin
   *
   * @param reader the reader
   * @return the window
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private static int readGenesAddress(RandomAccessFile reader) throws IOException {
    reader.seek(GENES_BYTES_OFFSET);

    return reader.readInt();
  }

  /**
   * Convert a genomic type to a int representation.
   * 
   * @param type
   * @return
   */
  public static int getType(GenomicType type) {
    switch (type) {
    case TRANSCRIPT:
      return GENOMIC_TYPE_TRANSCRIPT;
    case EXON:
      return GENOMIC_TYPE_EXON;
    default:
      return GENOMIC_TYPE_GENE;
    }
  }

  /**
   * Gets the bin address.
   *
   * @param reader the reader
   * @param bin    the bin
   * @return the bin address
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private static int readBinAddress(RandomAccessFile reader, int bin) throws IOException {
    reader.seek(HEADER_BYTES_OFFSET + bin * INT_BYTES);

    return reader.readInt();
  }

  public static int getStrand(Strand strand) {
    if (Strand.isSense(strand)) {
      return 0;
    } else {
      return 1;
    }
  }

  public static Strand getStrand(int strand) {
    if (strand == 0) {
      return Strand.SENSE;
    } else {
      return Strand.ANTISENSE;
    }
  }

}