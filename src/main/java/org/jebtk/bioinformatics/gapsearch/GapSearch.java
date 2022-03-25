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

import java.util.Collections;
import java.util.List;
import java.util.Map.Entry;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.collections.UniqueArrayList;

/**
 * Generic interface for quickly searching for features by genomic location.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public abstract class GapSearch<T> implements Iterable<Chromosome> {

  /**
   * Adds the feature.
   *
   * @param region  the region
   * @param feature the feature
   */
  public abstract void add(GenomicRegion region, T feature);

  /**
   * Size.
   *
   * @return the int
   */
  public abstract int size();

  /**
   * Gets the feature list.
   *
   * @return the feature list
   */
  public abstract List<T> getFeatures();

  /**
   * Gets values in the blocks spanning the region of interest.
   *
   * @param region the region
   * @return the feature list
   */
  public List<T> getValues(GenomicRegion region) {
    List<GappedSearchFeatures<T>> range = getFeatures(region);

    if (range.size() == 0) {
      return Collections.emptyList();
    }

    List<T> ret = new UniqueArrayList<T>();

    for (GappedSearchFeatures<T> features : range) {
      for (Entry<GenomicRegion, List<T>> r : features) {
        ret.addAll(r.getValue());
      }
    }

    return ret;
  }

  /**
   * Gets the feature list.
   *
   * @param chr the chr
   * @return the feature list
   */
  public abstract List<T> getFeatures(Chromosome chr);

  public abstract boolean contains(Chromosome chr);

  /**
   * Gets the features.
   *
   * @param region the region
   * @return the features
   */
  public List<GappedSearchFeatures<T>> getFeatures(GenomicRegion region) {
    if (region == null) {
      return Collections.emptyList();
    }

    return getFeatures(region.getChr(), region.getStart(), region.getEnd());
  }

  /**
   * Gets the features.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the features
   */
  public abstract List<GappedSearchFeatures<T>> getFeatures(Chromosome chr, int start, int end);

  /**
   * Returns the features closest to the region.
   * 
   * @param region
   * @return
   */
  public abstract List<T> getClosestFeatures(GenomicRegion region);

  public List<T> find(GenomicRegion region, int minBp) {
    return getOverlappingFeatures(region, minBp).toList();
  }

  public SearchResults<T> getOverlappingFeatures(GenomicRegion region, int minBp) {
    SearchResults<T> ret = new SearchResults<T>();

    getOverlappingFeatures(region, minBp, ret);

    return ret;
  }

  public void getOverlappingFeatures(GenomicRegion region, int minBp, SearchResults<T> ret) {
    List<GappedSearchFeatures<T>> allFeatures = getFeatures(region);

    if (allFeatures.size() == 0) {
      return;
    }

    for (GappedSearchFeatures<T> features : allFeatures) {
      for (Entry<GenomicRegion, List<T>> r : features) {
        GenomicRegion overlap = GenomicRegion.overlap(region, r.getKey());

        if (overlap != null && overlap.getLength() >= minBp) {
          ret.addAll(r.getKey(), r.getValue());
        }
      }
    }
  }

  /**
   * Checks for overlapping features.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return true, if successful
   */
  public boolean hasOverlappingFeatures(GenomicRegion region, int minBp) {
    List<GappedSearchFeatures<T>> allFeatures = getFeatures(region);

    if (allFeatures.size() == 0) {
      return false;
    }

    for (GappedSearchFeatures<T> features : allFeatures) {
      for (Entry<GenomicRegion, List<T>> r : features) {
        GenomicRegion overlap = GenomicRegion.overlap(region, r.getKey());

        if (minBp == -1 || overlap.getLength() >= minBp) {
          return true;
        }
      }
    }

    return false;
  }
}
