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

import java.io.DataInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jebtk.bioinformatics.gapsearch.BinaryGapSearch;
import org.jebtk.bioinformatics.gapsearch.BinarySearch;
import org.jebtk.bioinformatics.gapsearch.GappedSearchFeatures;
import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.core.io.FileUtils;

/**
 * Fast search of genome sequence Paths to get get actual genomic data. This
 * Path reads 4bit encoded genomes (i.e. 2 bases per byte).
 *
 * @author Antony Holmes
 *
 */
public class SequenceReader2Bit extends ChrSequenceReader {

  /** The m N map. */
  private BinaryGapSearch<GenomicRegion> mNMap = new BinarySearch<GenomicRegion>(); // new
                                                                                    // BinarySearch<GenomicRegion>();

  /** The m mask map. */
  private BinaryGapSearch<GenomicRegion> mMaskMap = new BinarySearch<GenomicRegion>(); // new
                                                                                       // BinarySearch<GenomicRegion>();

  /** The m offset map. */
  private Map<Chromosome, Integer> mOffsetMap = new HashMap<Chromosome, Integer>(100);

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public SequenceReader2Bit(Path directory) {
    super(directory);
  }

  @Override
  public String getName() {
    return "2bit";
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
    Chromosome chr = region.getChr();

    if (!mFileMap.containsKey(chr)) {
      Path file = mFile.resolve(chr + ".2bit.gz");

      mFileMap.put(chr, file);

      loadMaskData(genome, chr, file);
    }

    return new SequenceRegion(region, getSequence2Bit(mFileMap.get(chr), genome, chr, region.mStart, region.mEnd,
        mOffsetMap.get(chr), displayUpper, repeatMaskType));
  }

  /**
   * Load mask data.
   *
   * @param chr  the chr
   * @param file the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private void loadMaskData(Genome genome, Chromosome chr, Path file) throws IOException {
    // TODO Auto-generated method stub

    DataInputStream in = FileUtils.newDataInputStream(file);

    try {
      int check = in.readInt();
      int version = in.readInt();
      int offset = in.readInt();

      mOffsetMap.put(chr, offset);

      int nc = in.readInt();

      for (int i = 0; i < nc; ++i) {
        GenomicRegion region = new GenomicRegion(chr, in.readInt(), in.readInt());

        mNMap.add(region, region);
      }

      int mc = in.readInt();

      for (int i = 0; i < mc; ++i) {
        GenomicRegion region = new GenomicRegion(chr, in.readInt(), in.readInt());

        mMaskMap.add(region, region);
      }
    } finally {
      in.close();
    }

  }

  /**
   * Read a sequence from a file assuming each base is encoded in 2 bits (4 bases
   * per byte).
   *
   * @param Path           the Path
   * @param chr            the chr
   * @param start          the start
   * @param end            the end
   * @param offset         the offset
   * @param displayUpper   the display upper
   * @param repeatMaskType the repeat mask type
   * @return the sequence4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public Sequence getSequence2Bit(Path Path, Genome genome, Chromosome chr, int start, int end, int offset,
      boolean displayUpper, RepeatMaskType repeatMaskType) throws IOException {

    int s = start - 1;
    int e = end - 1;

    byte[] buf = getBytes2Bit(Path, s, e, offset);

    // how many characters to read
    int l = end - start + 1;

    char[] buffer = new char[l];

    // byte mask;
    int v = 0;

    // the offset to start reading from
    // int b = s; // % 4;
    int bi = 0; // b / 4;
    int block;

    for (int i = 0; i < l; ++i) {
      block = s % 4;

      // System.err.println("b " + b + " " + buf[bi] + " " + v + " " + bi);

      switch (block) {
      case 0:
        v = ((buf[bi] & 192) >> 6);
        break;
      case 1:
        v = ((buf[bi] & 48) >> 4);
        break;
      case 2:
        v = ((buf[bi] & 12) >> 2);
        break;
      default:
        v = (buf[bi] & 3);
        // We are at the end of a byte so the next read must skip to
        // the next byte in the array
        ++bi;
        break;
      }

      char c = toChar(v); // , repeatMaskType);

      GenomicRegion testRegion = new GenomicRegion(chr, s, s);

      //
      // Determine if N
      //

      if (overlaps(mNMap, testRegion)) {
        c = 'N';
      }

      //
      // Determine if masked
      //

      // System.err.println("test repeat region " + testRegion);

      if (repeatMaskType != RepeatMaskType.UPPERCASE) {
        if (overlaps(mMaskMap, testRegion)) {
          if (repeatMaskType == RepeatMaskType.LOWERCASE) {
            c = toLower(c);
          } else {
            // Mask to N
            c = 'N';
          }
        }
      }

      buffer[i] = c;

      ++s;
    }

    String dna = new String(buffer);

    if (displayUpper) {
      return Sequence.create(dna);
    } else {
      return Sequence.create(dna.toLowerCase());
    }
  }

  /**
   * Determine if region overlaps an N or mask region.
   *
   * @param binSearch  the bin search
   * @param testRegion the test region
   * @return true, if successful
   */
  private static boolean overlaps(BinaryGapSearch<GenomicRegion> binSearch, GenomicRegion testRegion) {
    List<GappedSearchFeatures<GenomicRegion>> allFeatures = binSearch.getFeatures(testRegion);

    for (GappedSearchFeatures<GenomicRegion> feature : allFeatures) {
      for (Entry<GenomicRegion, List<GenomicRegion>> r : feature) {
        if (GenomicRegion.overlaps(r.getKey(), testRegion)) {
          return true;
        }
      }
    }

    return false;
  }

  /**
   * Gets the bytes4 bit.
   *
   * @param file   the file
   * @param start  the start
   * @param end    the end
   * @param offset the offset
   * @return the bytes4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static byte[] getBytes2Bit(Path file, int start, int end, int offset) throws IOException {
    int sb = start / 4 + offset;
    int eb = end / 4 + offset;

    // System.err.println(sb + " " + eb);

    return getBytes(file, sb, eb);
  }
}
