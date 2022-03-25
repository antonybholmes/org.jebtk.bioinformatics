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
package org.jebtk.bioinformatics.conservation;

import java.io.IOException;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.FileSequenceReader;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.io.PathUtils;

/**
 * Stores 1 bit per genome base ergo this is only suitable for binary
 * representations of the genome e.g. base matches other species, true or false.
 *
 * @author Antony Holmes
 *
 */
public class AlignmentAssemblyFile1Bit extends ConservationAssembly {

  /**
   * The member directory.
   */
  protected Path mDirectory;

  /**
   * The member file map.
   */
  protected Map<Chromosome, Path> mFileMap = new HashMap<Chromosome, Path>();

  /**
   * Directory containing genome files. Each file must contain exactly one line
   * consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public AlignmentAssemblyFile1Bit(Path directory) {
    mDirectory = directory;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.conservation.ConservationAssembly#
   * getScores(edu.columbia.rdf.lib.bioinformatics.genome.GenomicRegion)
   */
  @Override
  public List<Double> getScores(GenomicRegion region) throws IOException {
    Chromosome chr = region.getChr();

    if (!mFileMap.containsKey(chr)) {
      mFileMap.put(chr, mDirectory.resolve(chr + ".1bit.hg19.mm10"));
    }

    return getScores(mFileMap.get(chr), region.getStart(), region.getEnd());
  }

  /**
   * Gets the scores.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return the scores
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Double> getScores(Path file, int start, int end) throws IOException {

    int s = start - 1;
    int e = end - 1;

    byte[] buf = getBytes1Bit(file, s, e);

    // how many characters to read
    int l = end - start + 1;

    List<Double> scores = new ArrayList<Double>(l);

    // byte mask;
    int v;

    // the offset to start reading from
    int b = s % 8;

    for (int i = 0; i < l; ++i) {
      int bi = b / 8;

      int offset = b % 8;

      switch (offset) {
      case 0:
        v = (buf[bi] & 128) >> 7;
        break;
      case 1:
        v = (buf[bi] & 64) >> 6;
        break;
      case 2:
        v = (buf[bi] & 32) >> 5;
        break;
      case 3:
        v = (buf[bi] & 16) >> 4;
        break;
      case 4:
        v = (buf[bi] & 8) >> 3;
        break;
      case 5:
        v = (buf[bi] & 4) >> 2;
        break;
      case 6:
        v = (buf[bi] & 2) >> 1;
        break;
      default:
        v = buf[bi] & 1;
        break;
      }

      scores.add((double) v);

      ++b;
    }

    return scores;
  }

  /**
   * Gets the bytes1 bit.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return the bytes1 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static byte[] getBytes1Bit(Path file, int start, int end) throws IOException {
    int sb = start / 8;
    int eb = end / 8;

    return FileSequenceReader.getBytes(file, sb, eb);
  }

  /**
   * The main method.
   *
   * @param args the arguments
   * @throws IOException    Signals that an I/O exception has occurred.
   * @throws ParseException the parse exception
   */
  public static void main(String[] args) throws IOException, ParseException {
    AlignmentAssemblyFile1Bit a = new AlignmentAssemblyFile1Bit(
        PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/alignment/human_mouse_hg19_mm10"));

    // System.err.println(a.getScores("mm10", "chr1:11872-12139"));

    // System.err.println(a.getScores("chr10:87575-87801") + " " +
    // Statistics.pNonZero(a.getScores("chr10:87575-87801")));
  }
}
