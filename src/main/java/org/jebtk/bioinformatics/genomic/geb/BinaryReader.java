package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.MMapReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

abstract class BinaryReader {

  protected static final Logger LOG = LoggerFactory.getLogger(BinaryReader.class);

  private final Path mDir;
  private MMapReader mReader;
  protected final Genome mGenome;
  protected final int mWindow;
  private final Chromosome mChr;

  protected final String mPrefix;

  public BinaryReader(Path dir, String prefix, Genome genome, int window) {
    this(dir, prefix, genome, null, window);
  }

  public BinaryReader(Path dir, String prefix, Genome genome, Chromosome chr, int window) {
    mDir = dir;
    mPrefix = prefix;
    mGenome = genome;
    mChr = chr;
    mWindow = window;
  }

  public MMapReader getReader() throws IOException {
    return getReader(mChr);
  }

  public MMapReader getReader(Chromosome chr) throws IOException {
    /*
     * if (chr != null) { if (mChr == null || !chr.equals(mChr)) { if (mReader !=
     * null) { mReader.close(); }
     * 
     * mReader = FileUtils.newMemMappedReader(mDir.resolve(getFileName(chr))); mChr
     * = chr; } } else { if (mReader == null) { mReader =
     * FileUtils.newMemMappedReader(mDir.resolve(getFileName(chr))); } }
     */

    if (mReader == null) {
      Path file = mDir.resolve(getFileName(chr));
      LOG.info("Creating reader {}...", file);
      mReader = FileUtils.newMemMappedReader(file);
      // mReader = FileUtils.newRandomAccess(file);
    }

    return mReader;
  }

  public void close() throws IOException {
    if (mReader != null) {
      mReader.close();
    }
  }

  protected abstract Path getFileName(Chromosome chr);

  protected MMapReader seek(long address) throws IOException {
    return getReader().seek(address);
  }

  protected long tell() throws IOException {
    return getReader().tell();
  }

  protected int read() throws IOException {
    return getReader().read();
  }

  protected int readInt() throws IOException {
    return getReader().readInt();
  }

  protected double readDouble() throws IOException {
    return getReader().readDouble();
  }

  /*
   * protected Chromosome readChr() throws IOException { int c =
   * getReader().read();
   * 
   * String chr;
   * 
   * if (c < 32) { chr = "chr" + c; } else { switch(c) { case 30: chr = "chrX";
   * break; case 31: chr = "chrY"; default: chr = "chrM"; break; } }
   * 
   * return Chromosome.newChr(chr, mGenome); }
   */

  /**
   * Return the check number of 42 to indicate file is being read correctly. If 42
   * is not returned, the file is corrupt.
   * 
   * @return
   * @throws IOException
   */
  public int readCheckNum() throws IOException {
    return seek(0).readInt();
  }
}
