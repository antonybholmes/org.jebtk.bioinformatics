package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.collections.ArrayUtils;
import org.jebtk.core.collections.UniqueArrayList;

public class BinReader extends BinaryReader {
  private static final int MIN_BYTES_OFFSET = GEBReader.WINDOW_BYTE_OFFSET + GEBReader.INT_BYTES;

  private static final int N_BYTES_OFFSET = MIN_BYTES_OFFSET + GEBReader.INT_BYTES;

  public static final int HEADER_BYTES_OFFSET = N_BYTES_OFFSET + GEBReader.INT_BYTES;

  private int mMinBin = -1;

  public BinReader(Path dir, String prefix, Genome genome, Chromosome chr, int window) throws IOException {
    super(dir, prefix, genome, chr, window);
  }

  // private static int convertAddress(int address) {
  // return HEADER_BYTES_OFFSET + address;
  // }

  // public void seek(int address) throws IOException {
  // getReader().seek(convertAddress(address));
  // }

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
  private int[] binAddressesFromBins(final int[] bins) throws IOException {
    int n = bins.length;

    int[] ret = new int[n];

    for (int i = 0; i < n; ++i) {
      ret[i] = readBinAddress(bins[i]);

      LOG.info("Reading bin from {} {} {}", i, bins[i], ret[i]);
    }

    return ret;
  }

  private int minBin() throws IOException {
    if (mMinBin == -1) {
      seek(MIN_BYTES_OFFSET);
      mMinBin = readInt();
    }

    return mMinBin;
  }

  /**
   * Gets the bin address.
   *
   * @param reader the reader
   * @param bin    the bin
   * @return the bin address
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int readBinAddress(int bin) throws IOException {
    if (mMinBin == -1) {
      seek(MIN_BYTES_OFFSET);
      seek(MIN_BYTES_OFFSET);
      mMinBin = readInt();
      System.err.println("min bin " + mMinBin);
    }

    seek(HEADER_BYTES_OFFSET + (bin - mMinBin) * GEBReader.INT_BYTES);

    return readInt();
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
  private List<Integer> elementAddressesFromBins(int[] binAddresses) throws IOException {

    // Gene address
    int ga;

    // How many genes are in the bin
    int size;

    List<Integer> ret = new UniqueArrayList<Integer>(binAddresses.length);

    for (int ba : binAddresses) {
      seek(ba);

      // 255
      // getReader().read();

      size = getReader().readInt();

      System.err.println("size up " + ba + " " + size);

      // Read how many addresses are in the bin and then extract them
      for (int j = 0; j < size; ++j) {
        ga = getReader().readInt();

        ret.add(ga);
      }
    }

    // Return the number of genes we found
    return ret;
  }

  public List<Integer> elementAddresses(GenomicRegion region) throws IOException {
    int sb = region.mStart / mWindow;
    int eb = region.mEnd / mWindow;

    int[] bins = ArrayUtils.array(sb, eb);

    int[] binAddresses = binAddressesFromBins(bins);

    System.err.println(Arrays.toString(bins) + " " + Arrays.toString(binAddresses));

    List<Integer> elementAddresses = elementAddressesFromBins(binAddresses);

    return elementAddresses;
  }

  @Override
  protected Path getFileName(Chromosome chr) {
    return getFileName(mPrefix, chr);
  }

  public static final Path getFileName(String prefix, Chromosome chr) {
    return GEBReader.getFileName("bins", prefix, chr);
  }
}
