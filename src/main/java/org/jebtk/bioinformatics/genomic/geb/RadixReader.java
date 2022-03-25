package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;

public class RadixReader extends BinaryReader {
  public static final int HEADER_BYTES_OFFSET = GEBReader.WINDOW_BYTE_OFFSET + GEBReader.INT_BYTES;

  public static final int RADIX_TREE_PREFIX_BYTES = 1 + GEBReader.INT_BYTES;

  public RadixReader(Path dir, String prefix, Genome genome, int window) throws IOException {
    super(dir, prefix, genome, window);
  }

  @Override
  protected Path getFileName(Chromosome chr) {
    return getFileName(mPrefix);
  }

  public List<Integer> elementAddresses(String id) throws IOException {
    return elementAddresses(id, false);
  }

  public List<Integer> elementAddresses(String id, boolean exact) throws IOException {
    List<Integer> ret = new ArrayList<Integer>();

    elementAddresses(id, exact, ret);

    return ret;
  }

  public void elementAddresses(String id, boolean exact, List<Integer> ret) throws IOException {

    char[] ca = id.toLowerCase().toCharArray();

    // Find the tree start
    seek(HEADER_BYTES_OFFSET);

    char leafc = 0;
    int address = 0;
    int n = 0;

    boolean found = false;

    // int m = 0;

    for (char c : ca) {

      // Number of children
      n = read(); // .readInt();

      // assume we won't find a match
      found = false;

      for (int i = 0; i < n; ++i) {
        leafc = (char) read(); // readByte();
        address = readInt();

        if (leafc == c) {
          // we did find a match so keep going
          found = true;
          seek(address);
          break;
        }
      }

      if (!found) {
        break;
      }
    }

    if (!found) {
      return;
    }

    // This means we kept finding a child matching the prefix so the seek
    // is at the beginning of a node either because we ran out of chars
    // or nodes. In this case we must skip over the child addresses and
    // just look at the addresses of the genes

    // skip past number of addresses and the addresses themselves
    // getReader().seek(address + N_BYTES + getReader().readInt() *
    // RADIX_TREE_PREFIX_BYTES);
    getReader().skipBytes(read() * RADIX_TREE_PREFIX_BYTES); // readInt()

    // Should be on a node that is hopefully matches our search term
    // Since we checked all the children, the seek is at the position of
    // the gene addressses so we can get them

    n = readInt();

    for (int i = 0; i < n; ++i) {
      ret.add(readInt());
    }

    if (exact) {
      return;
    }

    // Add the partial
    n = readInt();

    for (int i = 0; i < n; ++i) {
      ret.add(readInt());
    }
  }

  public static final Path getFileName(String prefix) {
    return GEBReader.getFileName("radix", prefix);
  }
}
