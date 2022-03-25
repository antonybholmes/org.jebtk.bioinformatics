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
import java.nio.file.Path;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceRegion;

/**
 * Fast search of genome sequence Paths to get get actual genomic data. This
 * Path reads 4bit encoded genomes (i.e. 2 bases per byte).
 *
 * @author Antony Holmes
 *
 */
public class SequenceReader4Bit extends ChrSequenceReader {

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public SequenceReader4Bit(Path directory) {
    super(directory);
  }

  @Override
  public String getName() {
    return "4bit";
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
      mFileMap.put(chr, mFile.resolve(chr + ".4bit"));
    }

    return new SequenceRegion(region,
        getSequence4Bit(mFileMap.get(chr), region.getStart(), region.getEnd(), displayUpper, repeatMaskType));
  }

  /**
   * Gets the sequence4 bit.
   *
   * @param Path           the Path
   * @param start          the start
   * @param end            the end
   * @param displayUpper   the display upper
   * @param repeatMaskType the repeat mask type
   * @return the sequence4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Sequence getSequence4Bit(Path Path, int start, int end, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {

    int s = start - 1;
    int e = end - 1;

    byte[] buf = getBytes4Bit(Path, s, e);

    // how many characters to read
    int l = end - start + 1;

    char[] buffer = new char[l];

    // byte mask;
    int v;

    // the offset to start reading from
    int b = s % 2;
    int bi = b / 2;

    for (int i = 0; i < l; ++i) {
      if (b % 2 == 0) {
        // The upper 4 bits as a mask
        v = (buf[bi] & 240) >> 4;
      } else {
        v = buf[bi] & 15;
        ++bi;
      }

      // v = buf[bi] & mask;

      // System.err.println("b " + b + " " + buf[bi] + " " + v + " " + bi);

      buffer[i] = toChar(v, repeatMaskType);

      ++b;
    }

    String dna = new String(buffer);

    if (displayUpper) {
      return Sequence.create(dna);
    } else {
      return Sequence.create(dna.toLowerCase());
    }
  }

  /**
   * Gets the bytes4 bit.
   *
   * @param Path  the Path
   * @param start the start
   * @param end   the end
   * @return the bytes4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static byte[] getBytes4Bit(Path Path, int start, int end) throws IOException {
    int sb = start / 2;
    int eb = end / 2;

    // System.err.println(sb + " " + eb);

    return getBytes(Path, sb, eb);
  }
}
