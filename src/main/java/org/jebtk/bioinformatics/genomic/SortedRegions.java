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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.CollectionUtils;

/**
 * Keeps together sorted regions on the same chromosome. Regions are auto-sorted
 * before search
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class SortedRegions<T extends GenomicRegion> implements Iterable<T> {
  /**
   * Keep the regions sorted.
   */
  private List<T> mRegions = new ArrayList<T>(100);

  /**
   * The member sorted.
   */
  private boolean mSorted = false;

  /**
   * Adds the.
   *
   * @param region the region
   */
  public void add(T region) {
    mRegions.add(region);

    mSorted = false;
  }

  /**
   * Sort regions by start position if they are not sorted.
   */
  private void autoSort() {
    if (mSorted) {
      return;
    }

    mRegions = GenomicRegion.sortSingleChr(mRegions);

    mSorted = true;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<T> iterator() {
    autoSort();

    return mRegions.iterator();
  }

  /**
   * Size.
   *
   * @return the int
   */
  public int size() {
    return mRegions.size();
  }

  /**
   * Given a range, returns elements within that range (inclusive of start and
   * end).
   *
   * @param start the start
   * @param end   the end
   * @return the list
   */
  public List<T> findRegionsInRange(int start, int end) {
    autoSort();

    int s = findStartIndex(start);

    int l = Mathematics.getLength(s, findEndIndex(end));

    // System.err.println("find reg l " + l);

    return CollectionUtils.subList(mRegions, s, l);
  }

  /**
   * Find start index.
   *
   * @param start the start
   * @return the int
   */
  private int findStartIndex(int start) {
    if (mRegions.size() == 0) {
      return -1;
    }

    if (mRegions.get(mRegions.size() - 1).getStart() < start) {
      return -1;
    }

    if (start < mRegions.get(0).getStart()) {
      return 0;
    }

    // since the indices are ordered by start find
    // the smallest start greater than this start

    int s = 0;
    int e = mRegions.size() - 1;
    int mid;

    while (e - s > 1) {
      mid = (s + e) / 2;

      if (mRegions.get(mid).getStart() < start) {
        s = mid;
      } else {
        e = mid;
      }
    }

    return e;
  }

  /**
   * Find end index.
   *
   * @param end the end
   * @return the int
   */
  private int findEndIndex(int end) {
    if (mRegions.size() == 0) {
      return -1;
    }

    if (end < mRegions.get(0).getStart()) {
      return -1;
    }

    if (mRegions.get(mRegions.size() - 1).getEnd() < end) {
      return mRegions.size() - 1;
    }

    // since the indices are ordered by start find
    // the smallest start greater than this start

    int s = 0;
    int e = mRegions.size() - 1;
    int mid;

    while (e - s > 1) {
      mid = (s + e) / 2;

      if (mRegions.get(mid).getStart() < end) {
        s = mid;
      } else {
        e = mid;
      }
    }

    return s;
  }
}
