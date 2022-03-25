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
 * Fast search of genome sequence files to get get actual genomic data. This
 * file reads 4bit encoded genomes (i.e. 2 bases per byte).
 *
 * @author Antony Holmes
 *
 */
public class ConservationAssemblyFile8Bit extends ConservationAssembly {

  /**
   * The member directory.
   */
  protected Path mDirectory;

  /**
   * The member file map.
   */
  protected Map<Chromosome, Path> mFileMap = new HashMap<Chromosome, Path>();

  /**
   * Directory containing genome files which must be of the form chr.n.txt. Each
   * file must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public ConservationAssemblyFile8Bit(Path directory) {
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
      mFileMap.put(chr, mDirectory.resolve(chr + ".8bit.phastcons"));
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

    byte[] buf = FileSequenceReader.getBytes(file, s, e);

    List<Double> scores = new ArrayList<Double>(buf.length);

    for (int b : buf) {
      scores.add((double) b);
    }

    return scores;
  }

  /**
   * The main method.
   *
   * @param args the arguments
   * @throws IOException    Signals that an I/O exception has occurred.
   * @throws ParseException the parse exception
   */
  public static void main(String[] args) throws IOException, ParseException {
    ConservationAssemblyFile8Bit a = new ConservationAssemblyFile8Bit(
        PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/phastcons/hg19"));

    // System.err.println(a.getScores("hg19", "chr1:90040-90050"));
  }
}
