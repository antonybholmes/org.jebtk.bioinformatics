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
import java.util.Collection;
import java.util.List;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;

public abstract class DirsSequenceReader extends SequenceReader {

  /** The m directory. */
  protected final List<Path> mDirs = new ArrayList<Path>();

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public DirsSequenceReader(Path dir, Path... dirs) {
    mDirs.add(dir);

    for (Path d : dirs) {
      mDirs.add(d);
    }
  }

  public DirsSequenceReader(Path dir, Collection<Path> dirs) {
    mDirs.add(dir);

    mDirs.addAll(dirs);
  }

  @Override
  public String getName() {
    return "dirs";
  }

  public Path getDir() {
    return mDirs.get(0);
  }

  /**
   * Return the directories to search for assembly files.
   * 
   * @return
   */
  public Iterable<Path> getDirs() {
    return mDirs;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genome.GenomeAssembly#getGenomes()
   */
  @Override
  public List<Genome> getGenomes() throws IOException {

    List<Genome> ret = new ArrayList<Genome>();

    for (Path dir : mDirs) {
      List<Path> subDirs = FileUtils.lsdir(dir);

      for (Path sd : subDirs) {
        ret.add(GenomeService.getInstance().guessGenome(PathUtils.getName(sd)));
      }
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
    return getReader(genome).getSequence(genome, region, displayUpper, repeatMaskType);
  }

  @Override
  public List<SequenceRegion> getSequences(Genome genome, Collection<GenomicRegion> regions, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {
    return getReader(genome).getSequences(genome, regions, displayUpper, repeatMaskType);
  }

  public abstract SequenceReader getReader(Genome genome);
}
