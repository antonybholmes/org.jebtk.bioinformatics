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
package org.jebtk.bioinformatics.genomic;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.TextUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Deals with functions related to chromosomes.
 *
 * @author Antony Holmes
 *
 */
public class GenomeService {
  /**
   * The Class ChromosomesLoader.
   */
  private static class GenomeLoader {

    /** The Constant INSTANCE. */
    private static final GenomeService INSTANCE = new GenomeService();
  }

  /**
   * Gets the single instance of GenomeService.
   *
   * @return single instance of GenomeService.
   */
  public static GenomeService getInstance() {
    return GenomeLoader.INSTANCE;
  }

  private static final Logger LOG = LoggerFactory.getLogger(GenomeService.class);

  private GenomeGuess mGenomeGuess = new GenomeGuess();

  private Map<String, Genome> mGenomeMap = new HashMap<String, Genome>();

  /**
   * Instantiates a new chromosomes.
   */
  private GenomeService() {
    mGenomeMap.put(Genome.HG19.toString(), Genome.HG19);
    mGenomeMap.put(Genome.HUMAN_HG19_REFSEQ.toString(), Genome.HUMAN_HG19_REFSEQ);
    mGenomeMap.put(Genome.GRCH37.toString(), Genome.GRCH37);
    mGenomeMap.put(Genome.HUMAN_HG19_GENCODE.toString(), Genome.HUMAN_HG19_GENCODE);
    mGenomeMap.put(Genome.GRCH38.toString(), Genome.GRCH38);
    mGenomeMap.put(Genome.MM10.toString(), Genome.MM10);
    mGenomeMap.put(Genome.FVBJ.toString(), Genome.FVBJ);
  }

  public void add(Genome genome) {
    mGenomeMap.put(genome.toString(), genome);
  }

  /**
   * Returns a given genome object.
   * 
   * @param name
   * @return
   */
  public Genome get(String name) {
    return get(name, TextUtils.NA, TextUtils.NA);
  }

  /**
   * Returns a given genome object.
   * 
   * @param name
   * @param assembly
   * @return
   */
  public Genome get(String name, String assembly) {
    return get(name, assembly, TextUtils.NA);
  }

  /**
   * Returns a given genome.
   * 
   * @param name
   * @param assembly
   * @param track
   * @return
   */
  public Genome get(String name, String assembly, String track) {
    String id = Genome.genomeId(name, assembly, track);

    if (!mGenomeMap.containsKey(id)) {
      mGenomeMap.put(id, new Genome(name, assembly, track));
    }

    return mGenomeMap.get(id);
  }

  /**
   * Set the object for guessing the genome from an id.
   * 
   * @param guess
   */
  public void setGenomeGuess(GenomeGuess guess) {
    mGenomeGuess = guess;
  }

  public Genome guessGenome(String name) {
    return mGenomeGuess.guess(name);
  }

  public Genome guessGenome(Path file) {
    return guessGenome(PathUtils.getName(file));
  }
}
