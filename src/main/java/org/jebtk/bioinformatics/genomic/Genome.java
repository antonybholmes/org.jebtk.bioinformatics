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
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.jebtk.core.AppService;
import org.jebtk.core.NameGetter;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.DefaultTreeMapCreator;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.TreeMapCreator;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

import com.fasterxml.jackson.annotation.JsonGetter;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;

/**
 * The enum Genome.
 */
 @JsonPropertyOrder({ "name", "assembly", "source" })
public class Genome implements NameGetter, Comparable<Genome> {

  public static final Genome NO_GENOME = new Genome(TextUtils.NA);

  private static final String GENCODE = "gencode";

  private static final String REFSEQ = "refseq";

  public static final String HUMAN = "Human";
  public static final String MOUSE = "Mouse";

  public static final Genome HUMAN_HG18_UCSC = new Genome(HUMAN, "hg18", REFSEQ);
  public static final Genome HUMAN_HG19_UCSC = new Genome(HUMAN, "hg19", REFSEQ);

  public static final Genome HUMAN_HG18_REFSEQ = new Genome(HUMAN, "hg18", REFSEQ);
  public static final Genome HUMAN_HG19_REFSEQ = new Genome(HUMAN, "hg19", REFSEQ);
  public static final Genome HUMAN_HG19_GENCODE = new Genome(HUMAN, "hg19", GENCODE);
  public static final Genome HUMAN_GRCH37_GENCODE = new Genome(HUMAN, "GRCh37", GENCODE);
  public static final Genome HUMAN_GRCH38_GENCODE = new Genome(HUMAN, "GRCh38", GENCODE);

  public static final Genome HG18 = new Genome(HUMAN, "hg18"); // HUMAN_HG18_REFSEQ;
  public static final Genome HG19 = new Genome(HUMAN, "hg19"); // HUMAN_HG19_REFSEQ;

  public static final Genome GRCH37 = new Genome(HUMAN, "GRCh37");;

  public static final Genome GRCH38 = new Genome(HUMAN, "GRCh38");

  /** The Constant MM10. */
  public static final Genome MM10 = new Genome(MOUSE, "mm10", REFSEQ);
  
  public static final Genome FVBJ = new Genome(MOUSE, "fvbj", REFSEQ);

  public static final Genome GRCM38 = new Genome(MOUSE, "GRCm38", GENCODE);

  public static final Path GENOME_HOME = AppService.RES_HOME.resolve("genomes");

  /** Local genome dir of app */
  public static final Path GENOME_DIR = AppService.RES_DIR.resolve("genomes");

  // private List<Path> mDirs = new ArrayList<Path>();

  // private ChromosomeDirs mChrs = null;

  public final String mName;

  public final String mAssembly;

  public final String mTrack;

  private String mS;

  public Genome(String name) {
    this(name, TextUtils.NA);
  }

  public Genome(String name, String assembly) {
    this(name, assembly, TextUtils.NA);
  }

  /**
   * 
   * @param name     e.g. 'human'
   * @param assembly e.g. 'grch38'
   * @param track    e.g. 'gencode'
   */
  public Genome(String name, String assembly, String track) {
    mName = TextUtils.sentenceCase(name);
    mAssembly = assembly;
    mTrack = track;
  }

  public Genome(Genome genome) {
    this(genome.mName, genome.mAssembly, genome.mTrack);
  }

  @Override
  @JsonGetter("name")
  public String getName() {
    return getGenome();
  }

  @JsonGetter("genome")
  public String getGenome() {
    return mName;
  }

  /**
   * Returns the assembly.
   * 
   * @return e.g. grch38
   */
  @JsonGetter("assembly")
  public String getAssembly() {
    return mAssembly;
  }

  /**
   * Returns the track.
   * 
   * @return e.g. 'gencode'.
   */
  @JsonGetter("track")
  public String getTrack() {
    return mTrack;
  }

  // public List<Path> getDirs() {
  // return mDirs;
  // }

  // public ChromosomeDirs chrs() {
  // if (mChrs == null) {
  // mChrs = new ChromosomeDirs(this, mDirs);
  // }
  //
  // return mChrs;
  // }

  /**
   * Invalidate the cache so it is rebuilt.
   */
  // public void cache() {
  // mChrs.cache();
  // mChrs = null;
  // }

  /*
   * public Chromosome chr(String chr) { return chrs().chr(chr); }
   * 
   * public Chromosome chr(int chr) { return chrs().chr(chr); }
   * 
   * public Chromosome chr(Chromosome chr) { return chr(chr.toString()); }
   * 
   * public Chromosome randChr() { return chrs().randChr(); }
   */

  @Override
  public String toString() {
    if (mS == null) {
      mS = genomeId(mName, mAssembly, mTrack);
    }

    return mS;
  }

  @Override
  public int compareTo(Genome g) {
    int ret = mName.toLowerCase().compareTo(g.mName.toLowerCase());

    if (ret == 0 && !mAssembly.equals(TextUtils.NA)) {
      ret = mAssembly.toLowerCase().compareTo(g.mAssembly.toLowerCase());

      if (ret == 0 && !mTrack.equals(TextUtils.NA)) {
        ret = mTrack.toLowerCase().compareTo(g.mTrack.toLowerCase());
      }
    }

    return ret;
  }

  @Override
  public boolean equals(Object o) {
    if (o instanceof Genome) {
      return compareTo((Genome) o) == 0;
    } else {
      return false;
    }
  }

  @Override
  public int hashCode() {
    return toString().hashCode();
  }

  /**
   * Change the track on a genome since genomes are immutable.
   * 
   * @param genome
   * @param track
   * @return
   */
  public static Genome changeTrack(Genome genome, String track) {
    return new Genome(genome.mName, genome.mAssembly, track);
  }

  /**
   * Create a genome from a colon separated id e.g. human:grch38:gencode
   * 
   * @param g
   * @return
   */
  public static Genome fromId(String g) {
    List<String> tokens = Splitter.onColon().text(g); // TextUtils.fastSplit(g,
    // TextUtils.COLON_DELIMITER);

    switch (tokens.size()) {
    case 2:
      return new Genome(tokens.get(0), tokens.get(1));
    case 3:
      return new Genome(tokens.get(0), tokens.get(1), tokens.get(2));
    default:
      return new Genome(tokens.get(0));
    }
  }

  public static String genomeId(String name, String assembly) {
    return genomeId(name, assembly, TextUtils.NA);
  }

  /**
   * Create a genome id from the genome name parts.
   * 
   * @param name
   * @param assembly
   * @param track
   * @return
   */
  public static String genomeId(String name, String assembly, String track) {
    List<String> ret = new ArrayList<String>(3);

    ret.add(name.toLowerCase());

    if (!assembly.equals(TextUtils.NA)) {
      ret.add(assembly.toLowerCase());

      if (!track.equals(TextUtils.NA)) {
        ret.add(track.toLowerCase());
      }
    }

    // return Join.onColon().values(ret).toString();

    return ret.stream().map(String::valueOf).collect(Collectors.joining(":"));
  }

  /**
   * Sort genomes by name, assembly, track to put in alphabetical order
   * 
   * @param genomes
   * @return
   * @throws IOException
   */
  public static IterMap<String, IterMap<String, IterMap<String, Genome>>> sort(Iterable<Genome> genomes) {
    IterMap<String, IterMap<String, IterMap<String, Genome>>> genomeMap = DefaultTreeMap
        .create(new DefaultTreeMapCreator<String, IterMap<String, Genome>>(new TreeMapCreator<String, Genome>()));

    for (Genome genome : genomes) {
      genomeMap.get(genome.getName()).get(genome.getAssembly()).put(genome.getTrack(), genome);
    }

    return genomeMap;
  }

  public static IterMap<String, IterMap<String, Genome>> sortByAssembly(Iterable<Genome> genomes) {
    IterMap<String, IterMap<String, Genome>> genomeMap = DefaultTreeMap.create(new TreeMapCreator<String, Genome>());

    for (Genome genome : genomes) {
      genomeMap.get(genome.getName()).put(genome.getAssembly(), genome);
    }

    return genomeMap;
  }

  /**
   * Return a genome containing only the assembly information. This is useful for
   * chr lookups where the information will not be tied to a track.
   * 
   * @param genome
   * @return
   */
  public static Genome assembly(Genome genome) {
    return assembly(genome.mName, genome.mAssembly);
  }

  public static Genome assembly(String name, String assembly) {
    return new Genome(name, assembly);
  }
}
