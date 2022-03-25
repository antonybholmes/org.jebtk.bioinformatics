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
package org.jebtk.bioinformatics.ext.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.core.io.Io;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Server for genome feature annotations.
 *
 * @author Antony Holmes
 *
 */
public class Cytobands {

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(Cytobands.class);

  /**
   * The member cytobands map.
   */
  private Map<Chromosome, List<Cytoband>> mCytobandsMap = new HashMap<Chromosome, List<Cytoband>>();

  /**
   * Instantiates a new cytobands.
   * 
   * @param genome
   *
   * @param reader the reader
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public Cytobands(Genome genome, BufferedReader reader) throws IOException {

    String line;

    // List<RefSeqGene> genes = new ArrayList<RefSeqGene>();

    // Map<Chromosome, Set<Integer>> positionMap =
    // new HashMap<Chromosome, Set<Integer>>();

    try {
      // skip header
      line = reader.readLine();

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        Cytoband cytoband = Cytoband.parse(genome, line);

        if (!mCytobandsMap.containsKey(cytoband.getChr())) {
          mCytobandsMap.put(cytoband.getChr(), new ArrayList<Cytoband>());
        }

        mCytobandsMap.get(cytoband.getChr()).add(cytoband);
      }
    } finally {
      reader.close();
    }
  }

  /**
   * Gets the cytobands.
   *
   * @param chr the chr
   * @return the cytobands
   */
  public List<Cytoband> getCytobands(Chromosome chr) {
    return mCytobandsMap.get(chr);
  }
}