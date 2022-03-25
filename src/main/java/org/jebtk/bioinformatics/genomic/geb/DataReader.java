package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;

public class DataReader extends BinaryReader {
  public static final int HEADER_BYTES_OFFSET = GEBReader.WINDOW_BYTE_OFFSET + GEBReader.INT_BYTES;

  /** The buffer. */
  private final byte[] mBuffer = new byte[256];

  public DataReader(Path dir, String prefix, Genome genome, int window) throws IOException {
    super(dir, prefix, genome, window);
  }

  public String readTag(int address) throws IOException {
    return readVarchar(address);
  }

  public String readVarchar(int address) throws IOException {
    seek(address);
    return readVarchar();
  }

  public double readDouble(int address) throws IOException {
    seek(address);
    return readDouble();
  }

  public int readInt(int address) throws IOException {
    seek(address);
    return readInt();
  }

  /**
   * Read a variable number of bytes to create a string.
   *
   * @param reader the reader
   * @return the string
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private String readVarchar() throws IOException {
    // First int tells us the length of the string
    int n = getReader().read();

    // Read n bytes into the buffer
    getReader().read(mBuffer, 0, n);

    // Create string from buffer
    String s = new String(mBuffer, 0, n, StandardCharsets.UTF_8);

    return s;
  }

  @Override
  protected Path getFileName(Chromosome chr) {
    return getFileName(mPrefix);
  }

  public static final Path getFileName(String prefix) {
    return GEBReader.getFileName("data", prefix);
  }
}
