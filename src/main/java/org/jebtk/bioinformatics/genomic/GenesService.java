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
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.io.PathUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Keep track of genes associated with genomes.
 *
 * @author Antony Holmes
 */
public class GenesService implements Iterable<Entry<Genome, GenesDB>> {

  /**
   * The Class GenesServiceLoader.
   */
  private static class GenesServiceLoader {

    /** The Constant INSTANCE. */
    private static final GenesService INSTANCE = new GenesService();
  }

  /**
   * Gets the single instance of SettingsService.
   *
   * @return single instance of SettingsService
   */
  public static GenesService getInstance() {
    return GenesServiceLoader.INSTANCE;
  }

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(GenesService.class);

  /**
   * The member symbol map.
   */
  private IterMap<Genome, GenesDB> mGenesMap = new IterTreeMap<Genome, GenesDB>();

  /**
   * Track dbs by genome
   */
  private IterMap<String, List<Genome>> mGenomeMap = DefaultTreeMap.create(new ArrayListCreator<Genome>());

  private Genome mCurrentDb;

  private GenomeDbGuess mDbGuess = new GenomeDbGuess();

  /**
   * Instantiates a new gene service.
   */
  private GenesService() {
    // do nothing
  }

  public GenomicElementsDB getGenes(Genome g) {
    System.err.println("Found genome in gene service " + g + " " + mGenesMap.containsKey(g) + " " + mGenesMap.keySet());

    if (mGenesMap.containsKey(g)) {
      return mGenesMap.get(g);
    } else {
      return GenomicElementsDB.EMPTY;
    }
  }

  public boolean contains(Genome g) {
    return mGenesMap.containsKey(g);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<Entry<Genome, GenesDB>> iterator() {
    return mGenesMap.iterator();
  }

  public void put(GenesDB genes) {
    for (Genome g : genes.getGenomes()) {
      put(g, genes);
    }
  }

  public void put(Genome g, GenesDB genes) {
    System.err.println("gene service " + g);
    
    mGenesMap.put(g, genes);

    mGenomeMap.get(g.getName()).add(g);
    mGenomeMap.get(g.getAssembly()).add(g);
    mGenomeMap.get(g.getTrack()).add(g);

    mCurrentDb = g;
  }

  public Genome getCurrentGenome() {
    return mCurrentDb;
  }

  /**
   * Gets the genomes.
   *
   * @return the genomes
   */
  public Iterable<Genome> getGenomes() {
    return mGenesMap.keySet();
  }

  public String guessDb(String name) {
    return mDbGuess.guess(name);
  }

  public String guessDb(Path file) {
    return guessDb(PathUtils.getName(file));
  }

  /**
   * Return a list of gene dbs for a given genome.
   * 
   * @param genome
   * @return
   */
  public Iterable<Genome> getGeneDbs(String genome) {
    return mGenomeMap.get(genome);
  }

  /**
   * Returns the first available gene database for a genome. This is a helper
   * method for when it is desirable to get gene metadata where position is not
   * important.
   * 
   * @param genome
   * @return
   */
  public Genome getFirstGeneDb(String genome) {
    if (mGenomeMap.get(genome).size() > 0) {
      return mGenomeMap.get(genome).get(0);
    } else {
      return null;
    }
  }
}