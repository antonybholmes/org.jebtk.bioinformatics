package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.bioinformatics.genomic.TagType;

public class ElementReader extends BinaryReader {
  public static final int N_BYTES_OFFSET = GEBReader.WINDOW_BYTE_OFFSET + GEBReader.INT_BYTES;

  public static final int HEADER_BYTES_OFFSET = N_BYTES_OFFSET + GEBReader.INT_BYTES;

  private DataReader mDataReader;

  public ElementReader(DataReader dataReader, Path dir, String prefix, Genome genome, int window) throws IOException {
    super(dir, prefix, genome, window);

    mDataReader = dataReader;
  }

  /**
   * Return the total number of primary elements (not including children).
   * 
   * @return
   * @throws IOException
   */
  private int getN() throws IOException {
    getReader().seek(N_BYTES_OFFSET);

    return readInt();
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
  public List<GenomicElement> readAll(GenomicType type) throws IOException {
    int n = getN();

    List<GenomicElement> ret = new ArrayList<GenomicElement>(n);

    for (int i = 0; i < n; ++i) {
      // Read all the genes
      readElement(type, ret);
    }

    return ret;
  }

  public List<GenomicElement> readElements(Collection<Integer> addresses, GenomicType type) throws IOException {
    List<GenomicElement> ret = new ArrayList<GenomicElement>(addresses.size());

    readElements(addresses, type, ret);

    return ret;
  }

  public void readElements(Collection<Integer> addresses, GenomicType type, List<GenomicElement> ret)
      throws IOException {

    for (int address : addresses) {
      // System.err.println("seeking " + address + " in " +
      // getFileName((Chromosome)null));
      seek(address);
      readElement(type, ret);
    }
  }

  private void readElement(GenomicType type, List<GenomicElement> ret) throws IOException {
    readElement(type, null, 0, ret);
  }

  private void readElement(GenomicType type, GenomicElement parent, int depth, List<GenomicElement> ret)
      throws IOException {
    GenomicElement element = readElement();

    boolean correctType = element.getType().equals(type);

    // If we haven't found what we are looking for, keep exploring the
    // children

    // Number of children
    int n = getReader().readShort();

    for (int i = 0; i < n; ++i) {
      // If type requested is not a gene, pass null to indicate that the
      // transcripts should not add themselves to the gene
      readElement(type, correctType ? element : null, depth + 1, ret);
    }

    if (parent != null) {
      parent.addChild(element);
    }

    if (correctType) {
      ret.add(element);
    }
  }

  /**
   * Read an element from file.
   * 
   * @param reader
   * @return
   * @throws IOException
   */
  private GenomicElement readElement() throws IOException {

    // Skip id (int)
    // readInt();

    // read block
    // read();

    int address = readInt();

    GenomicType t = GenomicType.parse(mDataReader.readVarchar(address));

    GenomicRegion l = readLocation();

    // System.err.println("type " + t + " " + l);

    Strand strand = readStrand();

    GenomicElement gene = new GenomicElement(t, l, strand);

    readProperties(gene);

    readTags(gene);

    return gene;
  }

  private GenomicRegion readLocation() throws IOException {

    Chromosome chr = readChr();

    int start = readInt();
    int end = readInt();

    return GenomicRegion.create(chr, start, end);
  }

  private Chromosome readChr() throws IOException {
    return Chromosome.newChr(mDataReader.readVarchar(readInt()));
  }

  private Strand readStrand() throws IOException {
    int strand = getReader().read();

    return getStrand(strand);
  }

  /**
   * Load tags associated with an entity.
   * 
   * @param reader
   * @param e
   * @return
   * @throws IOException
   */
  private int readTags(GenomicElement e) throws IOException {
    int n = getReader().read();

    int address;

    for (int i = 0; i < n; ++i) {
      address = readInt();

      e.addTag(mDataReader.readTag(address));
    }

    return n;
  }

  public int readProperties(GenomicElement e) throws IOException {
    int n = getReader().read();

    for (int i = 0; i < n; ++i) {
      readProperty(e);
    }

    return n;
  }

  private void readProperty(GenomicElement e) throws IOException {
    // Address of key name
    int nameAddress = readInt();

    // Address of value
    TagType propType = TagType.parse(getReader().read());
    int valueAddress = readInt();

    String name = mDataReader.readVarchar(nameAddress);

    switch (propType) {
    case DOUBLE:
      e.setProperty(name, mDataReader.readDouble(valueAddress));
      break;
    case INT:
      e.setProperty(name, mDataReader.readInt(valueAddress));
      break;
    default:
      e.setProperty(name, mDataReader.readVarchar(valueAddress));
      break;
    }
  }

  @Override
  protected Path getFileName(Chromosome chr) {
    return getFileName(mPrefix);
  }

  public static Strand getStrand(int strand) {
    if (strand == 0) {
      return Strand.SENSE;
    } else {
      return Strand.ANTISENSE;
    }
  }

  /**
   * Gets the file name.
   *
   * @param genome the genome
   * @param chr    the chr
   * @param window the window
   * @return the file name
   */
  public static final Path getFileName(String prefix) {
    return GEBReader.getFileName("elements", prefix);
  }
}
