package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.core.collections.UniqueArrayList;

public class BTreeReader extends BinaryReader {
  public static final int BINS_BYTES_OFFSET = GEBReader.WINDOW_BYTE_OFFSET + GEBReader.INT_BYTES;

  public static final int HEADER_BYTES_OFFSET = BINS_BYTES_OFFSET + GEBReader.INT_BYTES;

  public static final int BTREE_CHILD_ADDRESSES_BYTES = 2 * GEBReader.INT_BYTES;

  public BTreeReader(Path dir, String prefix, Genome genome, Chromosome chr, int window) throws IOException {
    super(dir, prefix, genome, chr, window);
  }

  /**
   * Gets the bin address.
   *
   * @param reader the reader
   * @param bin    the bin
   * @return the bin address
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int binAddressFromBin(int bin) throws IOException {
    seek(HEADER_BYTES_OFFSET);

    boolean search = true;

    while (search) {
      // node bin
      int b = readInt();

      // Address in bin array
      int ba = readInt();

      // Address of left child
      int la = readInt();

      // Address of right child
      int ra = readInt();

      // System.err.println("btree " + bin + " " + b + " " + ba + " " + la + " "
      // + ra);

      if (bin < b) {
        if (la > 0) {
          seek(la);
        } else {
          // We've no more children to explore so return nearest bin
          return ba; // -1;
        }
      } else if (bin > b) {
        if (ra > 0) {
          seek(ra);
        } else {
          // We've no more children to explore
          return ba; // -1;
        }
      } else {
        // Found the bin we are looking for, so read the address
        return ba;
      }
    }

    return -1;
  }

  private List<Integer> binAddressesFromTree(long address1, long address2) throws IOException {
    List<Integer> ret = new ArrayList<Integer>();

    seek(address1);

    System.err.print("add " + address1 + " " + address2);

    // keep going until we have got to the last address
    while (tell() <= address2) {
      ret.add(readInt());
    }

    return ret;

    /*
     * int n = treeBinAddresses.length;
     * 
     * int[] ret = new int[n];
     * 
     * for (int i = 0; i < n; ++i) { seek(treeBinAddresses[i]); ret[i] = readInt();
     * 
     * if (i==0 || ret[i] == 1) { System.err.println("alert " + i + " " +ret[i] +
     * " " + treeBinAddresses[i]); } }
     * 
     * return ret;
     */
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
  private List<Integer> elementAddressesFromBins(List<Integer> binAddresses) throws IOException {

    // Gene address
    int ga;

    // How many genes are in the bin
    int n;

    List<Integer> ret = new UniqueArrayList<Integer>(binAddresses.size());

    for (int ba : binAddresses) {

      seek(ba);

      n = readInt();

      // Read how many addresses are in the bin and then extract them
      for (int i = 0; i < n; ++i) {
        ga = readInt();
        ret.add(ga);
      }
    }

    // Return the number of genes we found
    return ret;
  }

  public List<Integer> elementAddresses(Chromosome chr, int start, int end) throws IOException {
    int sb = start / mWindow;
    int eb = end / mWindow;

    int b1 = binAddressFromBin(sb);
    int b2 = binAddressFromBin(eb);

    // int[] treeBinAddresses;

    if (b1 < 0 && b2 < 0) {
      return Collections.emptyList();
    }

    if (b1 < 0) {
      b1 = b2;
      // treeBinAddresses = new int[]{b2};
    } else if (b2 < 0) {
      // treeBinAddresses = new int[]{b1};
      b2 = b1;
    } else {
      // treeBinAddresses = ArrayUtils.array(b1, b2, 4);
    }

    // System.err.println(b1 + " " + b2); // + " " +
    // Arrays.toString(treeBinAddresses));

    List<Integer> binAddresses = binAddressesFromTree(b1, b2); // treeBinAddresses);

    List<Integer> elementAddresses = elementAddressesFromBins(binAddresses);

    return elementAddresses;
  }

  @Override
  protected Path getFileName(Chromosome chr) {
    return getFileName(mPrefix, chr);
  }

  public static final Path getFileName(String prefix, Chromosome chr) {
    return GEBReader.getFileName("btree", prefix, chr);
  }
}
