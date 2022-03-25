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

import org.jebtk.bioinformatics.gapsearch.BinaryGapSearch;
import org.jebtk.bioinformatics.gapsearch.FixedGapSearch;
import org.jebtk.bioinformatics.gapsearch.GapSearch;
import org.jebtk.core.model.ListModel;

/**
 * Keeps a sorted list of regions.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class GenomicRegions<T extends GenomicElement> extends ListModel<T> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * Gets the fixed gap search.
   *
   * @param <X>     the generic type
   * @param regions the regions
   * @return the fixed gap search
   */
  public static <X extends GenomicElement> GapSearch<X> getFixedGapSearch(Iterable<X> regions) {
    GapSearch<X> search = new FixedGapSearch<X>();

    for (X region : regions) {
      search.add(region, region);
    }

    return search;
  }

  /**
   * Gets the binary search.
   *
   * @param <X>     the generic type
   * @param regions the regions
   * @return the binary search
   */
  public static <X extends GenomicElement> BinaryGapSearch<X> getBinarySearch(Iterable<X> regions) {
    BinaryGapSearch<X> search = new BinaryGapSearch<X>();

    for (X region : regions) {
      search.add(region, region);
    }

    return search;
  }

  /*
   * protected Map<Chromosome, SortedRegions<T>> mChrRegions = new
   * TreeMap<Chromosome, SortedRegions<T>>();
   * 
   * protected List<T> mRegions = new ArrayList<T>();
   * 
   * public void add(T region) { mRegions.add(region);
   * getRegionsByChr(region.getChr()).add(region);
   * 
   * fireChanged(); }
   * 
   * public SortedRegions<T> getRegionsByChr(Chromosome chromosome) { if
   * (!mChrRegions.containsKey(chromosome)) { mChrRegions.put(chromosome, new
   * SortedRegions<T>()); }
   * 
   * return mChrRegions.get(chromosome); }
   * 
   * public T get(int i) { return mRegions.get(i); }
   * 
   * public int size() { return mRegions.size(); }
   * 
   * @Override public Iterator<Chromosome> iterator() { return
   * mChrRegions.keySet().iterator(); }
   */
}
