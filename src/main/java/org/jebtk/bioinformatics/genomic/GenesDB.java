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
import java.util.ArrayList;
import java.util.List;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public abstract class GenesDB extends GenomicElementsDB {

  private static final long serialVersionUID = 1L;

  /**
   * Return the RefSeq ids used to index these genes.
   *
   * @return the ref seq ids
   */
  public Iterable<String> getRefSeqIds() {
    return getIds(Gene.REFSEQ_ID);
  }

  public Iterable<String> getNames() throws IOException {
    return getIds(Gene.GENE_NAME);
  }

  /**
   * Filter a list of gene elements to match a given type.
   * 
   * @param features
   * @param type
   * @return
   * @throws IOException
   */
  public List<GenomicElement> filterByType(List<GenomicElement> features, GenomicType type) {
    List<GenomicElement> ret = new ArrayList<GenomicElement>(features.size());

    for (GenomicElement feature : features) {
      if (feature.mType.equals(type)) {
        ret.add(feature);
      }
    }

    return ret;
  }

}