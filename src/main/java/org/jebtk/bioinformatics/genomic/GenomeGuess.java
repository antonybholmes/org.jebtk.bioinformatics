/**
 * Copyright 2017 Antony Holmes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.jebtk.bioinformatics.genomic;

/**
 * The Class ChrSpeciesGuess guess the species from a name.
 */
public class GenomeGuess {

  /**
   * Guess.
   *
   * @param id the id
   * @return the string
   */
  public Genome guess(String id) {
    String lid = id.toLowerCase();
    
    if (lid.contains("grch38")) {
      return Genome.GRCH38;
    } else if (lid.contains("hg18")) {
      return Genome.HG18;
    } else if (lid.contains("hg19")) {
      return Genome.HG19;
    } else if (lid.contains("mm10")) {
      return Genome.MM10;
    } else if (lid.contains("grcm38")) {
      return Genome.GRCM38;
    } else if (lid.contains("human")) {
      return Genome.GRCH38;
    } else if (lid.contains("mm10")) {
      return Genome.MM10;
    } else if (lid.contains("mouse")) {
      return Genome.MM10;
    } else if (lid.contains("fvbj")) {
      return Genome.FVBJ;
    } else {
      return Genome.NO_GENOME;
    }
  }
}
