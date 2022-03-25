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
package org.jebtk.bioinformatics.gapsearch;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;

/**
 * Uses a binary search to identify the closest points to a genomic coordinate
 * and returns all objects in the range this spans.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class BinarySearch<T> extends BinaryGapSearch<T> {

  /**
   * Instantiates a new binary search.
   */
  public BinarySearch() {
    super(1);
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#addFeature(edu.
   * columbia.rdf.lib.bioinformatics.genome.Chromosome, int, int,
   * java.lang.Object)
   */
  @Override
  public void add(GenomicRegion region, T feature) {
    Chromosome chr = region.getChr();
    int start = region.getStart();
    int end = region.getEnd();

    if (!mFeatures.get(chr).containsKey(start)) {
      mFeatures.get(chr).put(start, new GappedSearchFeatures<T>(start));
    }

    mFeatures.get(chr).get(start).add(region, feature);

    if (!mFeatures.get(chr).containsKey(end)) {
      mFeatures.get(chr).put(end, new GappedSearchFeatures<T>(end));
    }

    mFeatures.get(chr).get(end).add(region, feature);

    ++mSize;

    // Indicate that indexes will need to be rebuilt before searching
    mAutoSorted = false;
  }
}
