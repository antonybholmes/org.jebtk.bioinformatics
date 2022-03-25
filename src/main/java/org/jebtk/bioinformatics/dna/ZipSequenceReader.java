/**
 * Copyright (C) 2016, Antony Holmes
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  3. Neither the name of copyright holder nor the names of its contributors 
 *     may be used to endorse or promote products derived from this software 
 *     without specific prior written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 */
package org.jebtk.bioinformatics.dna;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Encodes DNA in a 2 bit file representing ACGT. All other characters such as N
 * map to A. Bases are encoded in two bits, so 4 bases per byte. A = 0, C = 1, G
 * = 2, T = 3. Files can be accompanied by a corresponding n. Data is loaded in
 * memory to speed it up.
 * 
 *
 * @author Antony Holmes
 *
 */
public class ZipSequenceReader extends DNASequenceReader {

  private static interface Process1Bit {
    /**
     * Process a 1 bit mask.
     * 
     * @param i       The current index
     * @param v       The 1bit value from the corresponding 1bit buffer at index i.
     * @param charBuf The char buffer to be updated.
     */
    public void process(int i, int v, char[] charBuf);
  }

  private static final Process1Bit N_PROCESSOR = new Process1Bit() {
    @Override
    public void process(int i, int v, char[] charBuf) {
      if (v == 1) {
        charBuf[i] = 'N';
      }
    }
  };

  private static final Process1Bit LOWER_PROCESSOR = new Process1Bit() {
    @Override
    public void process(int i, int v, char[] charBuf) {
      if (v == 1) {
        charBuf[i] = toLower(charBuf[i]);
      }
    }
  };

  public static final Logger LOG = LoggerFactory.getLogger(ZipSequenceReader.class);

  // Use fixed size arrays to cache chromosome features. Arrays are set
  // to be larger than the amount of data that will be cached.

  // private char[] mDnaBuf = new char[300000000];

  private final byte[] mChrBuf = new byte[100000000];
  private final byte[] mChrMaskBuf = new byte[50000000];
  private final byte[] mChrNBuf = new byte[50000000];

  /**
   * Store read bytes. We assume fewer than 2 million bases will be read at once.
   */
  private final byte[] mBuf = new byte[500000];

  private static final int MAX_SIZE_BP = 1000000;

  private final char[] mCharBuf = new char[MAX_SIZE_BP];

  private Chromosome mChr;

  /** How many 'N' were read */
  private int mNL = 0;
  private int mML = 0;

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public ZipSequenceReader(Path zip) {
    super(zip);
  }

  @Override
  public String getName() {
    return "zip";
  }

  @Override
  public List<SequenceRegion> getSequences(Genome genome, Collection<GenomicRegion> regions, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {

    IterMap<Chromosome, List<GenomicRegion>> chrMap = DefaultTreeMap.create(new ArrayListCreator<GenomicRegion>());

    for (GenomicRegion region : regions) {
      chrMap.get(region.getChr()).add(region);
    }

    Map<GenomicRegion, SequenceRegion> mSeqMap = new HashMap<GenomicRegion, SequenceRegion>();

    for (Entry<Chromosome, List<GenomicRegion>> f : chrMap) {
      for (GenomicRegion region : f.getValue()) {
        SequenceRegion sequence = getSequence(genome, region, displayUpper, repeatMaskType);

        mSeqMap.put(region, sequence);
      }
    }

    List<SequenceRegion> ret = new ArrayList<SequenceRegion>(regions.size());

    for (GenomicRegion region : regions) {
      ret.add(mSeqMap.get(region));
    }

    return ret;
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.genome.GenomeAssembly#getSequence(edu.
   * columbia.rdf.lib.bioinformatics.genome.GenomicRegion, boolean,
   * edu.columbia.rdf.lib.bioinformatics.genome.RepeatMaskType)
   */
  @Override
  public final SequenceRegion getSequence(Genome genome, GenomicRegion region, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {
    return new SequenceRegion(region, getSequence2Bit(region, displayUpper, repeatMaskType));
  }

  /**
   * Gets the sequence4 bit.
   *
   * @param file           the Path
   * @param chr            the chr
   * @param start          the start
   * @param end            the end
   * @param displayUpper   the display upper
   * @param repeatMaskType the repeat mask type
   * @return the sequence4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public Sequence getSequence2Bit(GenomicRegion region, boolean displayUpper, RepeatMaskType repeatMaskType)
      throws IOException {

    Chromosome chr = region.getChr();
    int start = region.getStart();
    int end = region.getEnd();

    // Cache a chromosome in memory for speed. For optimum results,
    // sort batches of coordinates by chromosome so that they are
    // processed in memory as much as possible with fewer cache collisions.
    // We cache one chr at a time as a trade off between memory and speed.

    if (mChr == null || !chr.equals(mChr)) {
      LOG.info("Caching chromosome {}...", chr);

      readBases(mFile, chr, ".dna.2bit", mChrBuf);
      mNL = readBases(mFile, chr, ".n.1bit", mChrNBuf);
      mML = readBases(mFile, chr, ".mask.1bit", mChrMaskBuf);

      mChr = chr;
    }

    int s = start - 1;
    int e = end - 1;

    getBytes2Bit(mChrBuf, s, e, mBuf);

    // how many characters to read
    int l = Math.min(MAX_SIZE_BP, end - start + 1);

    // Convert bytes to chars
    process2bit(s, l, mBuf, mCharBuf);

    //
    // Deal with undefined bases
    //

    int n = getN(chr, start, end, mBuf);

    if (n > 0) {
      process1bit(s, l, mBuf, N_PROCESSOR, mCharBuf);
    }

    //
    // Repeat Mask
    //

    if (repeatMaskType != RepeatMaskType.UPPERCASE) {
      // If the repeat mask is uppercase, we do nothing. Only when
      // set to either N or lowercase is it worth checking the mask.

      n = getMask(chr, start, end, mBuf);

      if (n > 0) {
        if (repeatMaskType == RepeatMaskType.N) {
          // If mask set, change to 'N'
          process1bit(s, l, mBuf, N_PROCESSOR, mCharBuf);
        } else {
          // If mask set, change letter to lowercase
          process1bit(s, l, mBuf, LOWER_PROCESSOR, mCharBuf);

        }
      }
    }

    //
    // Finalize
    //

    // If not uppercase, convert to lowercase
    if (!displayUpper) {
      toLower(mCharBuf, l);
    }

    return Sequence.create(region.getLocation(), new String(mCharBuf, 0, l));
  }

  /**
   * Base data is encoded in bytes so that more than one base can be represented
   * by a byte.
   * 
   * @param chr
   * @param fileMap
   * @param buf
   * @return
   * @throws IOException
   */
  private static int readBases(Path zip, Chromosome chr, String ext, byte[] buf) throws IOException {
    int n = -1;

    String file = chr + ext;

    ZipFile zipFile = FileUtils.newZipFile(zip);

    ZipEntry zipEntry = zipFile.getEntry(file);

    if (zipEntry != null) {
      InputStream in = FileUtils.newBufferedInputStream(zipFile, zipEntry);

      try {
        n = in.read(buf);
      } finally {
        in.close();
      }
    }

    return n;
  }

  /**
   * Convert byte to DNA bases where the byte represents 4 bases.
   * 
   * @param s
   * @param l
   * @param buf
   * @param charBuf
   */
  private static void process2bit(int s, int l, final byte[] buf, char[] charBuf) {
    // byte mask;
    int v = 0;

    // the offset to start reading from
    int b = s; // % 4;
    int bi = 0; // b / 4;
    int block;

    for (int i = 0; i < l; ++i) {
      block = b % 4;

      switch (block) {
      case 0:
        v = (buf[bi] >> 6);
        break;
      case 1:
        v = (buf[bi] >> 4);
        break;
      case 2:
        v = (buf[bi] >> 2);
        break;
      default:
        v = buf[bi];
        // We are at the end of a byte so the next read must skip to
        // the next byte in the array
        ++bi;
        break;
      }

      // AND with 3 to get the lowest 2 bits
      v &= 3;

      charBuf[i] = toChar(v);

      ++b;
    }
  }

  /**
   * Iterate over 1 bit array, updating the char array.
   * 
   * @param s         Start
   * @param l         Length
   * @param buf       1 bit buffer
   * @param charBuf   char buffer
   * @param processor
   */
  private static void process1bit(int s, int l, byte[] buf, Process1Bit processor, char[] charBuf) {
    int bi = 0;
    int b = s;
    int v;
    int block;

    for (int i = 0; i < l; ++i) {
      block = b % 8;

      switch (block) {
      case 0:
        v = (buf[bi] >> 7);
        break;
      case 1:
        v = (buf[bi] >> 6);
        break;
      case 2:
        v = (buf[bi] >> 5);
        break;
      case 3:
        v = (buf[bi] >> 4);
        break;
      case 4:
        v = (buf[bi] >> 3);
        break;
      case 5:
        v = (buf[bi] >> 2);
        break;
      case 6:
        v = (buf[bi] >> 1);
        break;
      default:
        v = buf[bi];
        // We are at the end of a byte so the next read must skip to
        // the next byte in the array
        ++bi;
        break;
      }

      // We are only interested in the first bit
      v &= 1;

      // if (v == 1) {
      // charBuf[i] = 'N';
      // }

      processor.process(i, v, charBuf);

      ++b;
    }
  }

  /**
   * Returns the number of Ns in a range.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return
   * @return the n The number of bytes read or -1 if the N mask file does not
   *         exist.
   * 
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public int getN(Chromosome chr, int start, int end, byte[] ret) throws IOException {
    if (mNL == 0) {
      return -1;
    }

    int s = start - 1;
    int e = end - 1;

    return getBytes1Bit(mChrNBuf, s, e, ret);
  }

  /**
   * Returns the repeat mask for a range.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return
   * @return the mask
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private int getMask(Chromosome chr, int start, int end, byte[] ret) throws IOException {
    if (mML == 0) {
      return -1;
    }

    int s = start - 1;
    int e = end - 1;

    return getBytes1Bit(mChrMaskBuf, s, e, ret);
  }

  /**
   * Gets the bytes4 bit.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return
   * @return the bytes4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static int getBytes2Bit(final byte[] buf, int start, int end, byte[] ret) throws IOException {
    int sb = start / 4;
    int eb = end / 4;

    // System.err.println(sb + " " + eb);

    return getBytes(buf, sb, eb, ret);
  }

  /**
   * Gets the bytes 1 bit.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return
   * @return the bytes 1 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static int getBytes1Bit(final byte[] buf, int start, int end, byte[] ret) throws IOException {
    int sb = start / 8;
    int eb = end / 8;

    // System.err.println(sb + " " + eb);

    return getBytes(buf, sb, eb, ret);
  }

  public static int getBytes(final byte[] buf, int start, int end, byte[] ret) throws IOException {

    int l = end - start + 1;

    System.arraycopy(buf, start, ret, 0, l);

    return l;
  }
}
