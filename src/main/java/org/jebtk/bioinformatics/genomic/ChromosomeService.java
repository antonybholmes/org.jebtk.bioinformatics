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

import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;

import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Deals with functions related to chromosomes.
 *
 * @author Antony Holmes
 *
 */
public class ChromosomeService extends GenomeDirs implements Iterable<Genome> {
  /**
   * The Class ChromosomesLoader.
   */
  private static class GenomeLoader {

    /** The Constant INSTANCE. */
    private static final ChromosomeService INSTANCE = new ChromosomeService();
  }

  /**
   * Gets the single instance of GenomeService.
   *
   * @return single instance of GenomeService.
   */
  public static ChromosomeService getInstance() {
    return GenomeLoader.INSTANCE;
  }

  private static final Logger LOG = LoggerFactory.getLogger(ChromosomeService.class);

  // private static final String EXT2 = "genome.txt.gz";

  // private IterMap<String, IterMap<String, IterMap<String, Genome>>>
  // mGenomeCacheMap =
  // DefaultHashMap.create(new DefaultHashMapCreator<String, IterMap<String,
  // Genome>>(new HashMapCreator<String, Genome>()));

  // private IterMap<String, Genome> mGenomeMap = new IterTreeMap<String,
  // Genome>();

  private IterMap<Genome, ChromosomeReader> mChrsMap = new IterTreeMap<Genome, ChromosomeReader>();

  private boolean mAutoLoad = true;

  /**
   * Instantiates a new chromosomes.
   */
  public ChromosomeService() {
    // Do nothing

    this(Genome.GENOME_DIR); // Genome.GENOME_HOME, Genome.GENOME_DIR);
  }

  public ChromosomeService(Path... dirs) {
    super(dirs);
  }

  /**
   * Get a genome. Service will attempt to auto discover genome data and populate
   * the genome object. If no data exists, an empty genome will be created.
   * 
   * @param genome
   * @return
   */
  /*
   * public Genome genome(String genome, String db) { try { autoLoad(); } catch
   * (IOException e) { e.printStackTrace(); }
   * 
   * // String fg = formatKey(genome);
   * 
   * if (!mGenomeMap.containsKey(genome)) { // If the genome does not exist,
   * create one mGenomeMap.put(genome, new Genome(genome, db, mDirs)); }
   * 
   * return mGenomeMap.get(genome); }
   * 
   * public Genome genome(Genome genome) { return genome(genome.getName(),
   * genome.getBuild()); }
   * 
   * public Genome genome(String genome) { return genome(guessGenome(genome)); }
   */

  private void autoLoad() throws IOException {
    if (mAutoLoad) {
      LOG.info("Checking {} for genome info.", mDirs);

      mAutoLoad = false;

      for (Path dir : mDirs) {
        if (FileUtils.isDirectory(dir)) {
          LOG.info("Checking {} for genome info.", dir);

          // String db = PathUtils.getName(dir);

          for (Path genomeDir : FileUtils.lsdir(dir)) {
            System.err.println(genomeDir);

            if (FileUtils.isDirectory(genomeDir)) {
              String g = PathUtils.getName(genomeDir); // TextUtils.sentenceCase(PathUtils.getName(genomeDir));

              System.err.println(g);

              for (Path assemblyDir : FileUtils.lsdir(genomeDir)) {
                if (FileUtils.isDirectory(assemblyDir)) {
                  String assembly = PathUtils.getName(assemblyDir);

                  // if (!mGenomeMap.containsKey(db)) {
                  LOG.info("Discovered genome {} in {}.", assembly, assemblyDir);

                  Genome genome = GenomeService.getInstance().get(g, assembly);

                  // mGenomeMap.put(db, genome);

                  mChrsMap.put(genome, new ChromosomeDirs(genome, assemblyDir));
                  // }
                }
              }
            }
          }
        }
      }
    }
  }

  /**
   * Invalidate the cache so it will be rebuilt.
   */
  public void cache() {
    mAutoLoad = true;
  }

  private ChromosomeReader autoLoad(Genome genome) {
    genome = Genome.assembly(genome);

    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    if (!mChrsMap.containsKey(genome)) {
      mChrsMap.put(genome, new ChromosomeDirs(genome));
    }

    return mChrsMap.get(genome);
  }

  public Chromosome chr(String chr) {
    return chr(Genome.NO_GENOME, chr);
  }

  /**
   * Returns the chr from a given genome reference.
   * 
   * @param genome The genome, e.g. 'Human'
   * @param chr    The chromosome, e.g. 'chr1'
   * 
   * @return The chromosome object from the desired genome.
   */
  public Chromosome chr(Genome genome, String chr) {
    return autoLoad(genome).chr(chr);
  }

  public int size(Genome genome, String chr) {
    return autoLoad(genome).size(chr(genome, chr));
  }

  public int size(Genome genome, Chromosome chr) {
    return autoLoad(genome).size(chr);
  }

  public Chromosome hg19(String chr) {
    return chr(Genome.HG19, chr);
  }

  public Chromosome guessChr(Path file, String chr) {
    return chr(GenomeService.getInstance().guessGenome(file), chr);
  }

  public Chromosome randChr(String genome, String db) {
    return randChr(genome, db);
  }

  public Chromosome randChr(String genome) {
    return randChr(GenomeService.getInstance().guessGenome(genome));
  }

  public Chromosome randChr(Genome genome) {
    return mChrsMap.get(genome).randChr();
  }

  @Override
  public Iterator<Genome> iterator() {
    return mChrsMap.keySet().iterator();
  }
}
