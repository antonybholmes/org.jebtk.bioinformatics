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
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomeSequenceReader;
import org.jebtk.bioinformatics.genomic.GenomeService;
import org.jebtk.bioinformatics.genomic.SequenceReader;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;

/**
 * Encodes DNA in a 2 bit file representing ACGT. All other characters such as N
 * map to A. Bases are encoded in two bits, so 4 bases per byte. A = 0, C = 1, G
 * = 2, T = 3. Files can be accompanied by a corresponding n
 * 
 *
 * @author Antony Holmes
 *
 */
public class DirZipSequenceReader extends GenomeSequenceReader {

  public DirZipSequenceReader() {
    this(Genome.GENOME_HOME, Genome.GENOME_DIR);
  }

  public DirZipSequenceReader(Path dir, Path... dirs) {
    super(dir, dirs);
  }

  @Override
  public String getName() {
    return "zip-dir";
  }

  @Override
  public List<Genome> getGenomes() throws IOException {

    List<Genome> ret = new ArrayList<Genome>();

    for (Path dir : getDirs()) {
      List<Path> files = FileUtils.endsWith(dir, "dna.zip");

      for (Path file : files) {
        ret.add(GenomeService.getInstance().guessGenome(PathUtils.namePrefix(file)));
      }
    }

    return ret;
  }

  @Override
  public SequenceReader getReader(Genome genome) {
    if (!mGenomeMap.containsKey(genome)) {
      // Path d = dir.resolve(genome);

      for (Path dir : mDirs) {
        if (FileUtils.isDirectory(dir)) {
          Path d = dir.resolve(genome.getName());

          if (FileUtils.isDirectory(d)) {
            Path d2 = d.resolve(genome.getAssembly());

            if (FileUtils.isDirectory(d2)) {
              Path zip = d2.resolve(genome + ".dna.zip");

              if (FileUtils.isFile(zip)) {
                mGenomeMap.put(genome, new ZipSequenceReader(zip));
              }
            }
          }
        }
      }
    }

    return mGenomeMap.get(genome);
  }
}
