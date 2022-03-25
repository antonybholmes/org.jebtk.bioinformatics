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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.IterHashMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.UniqueArrayList;

/**
 * An extension of the fixed gap search that can perform a binary search on to
 * look for the closest items to a genomic location. This
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class BinaryGapSearch<T> extends FixedGapSearch<T> {

  /**
   * The member auto sorted.
   */
  protected boolean mAutoSorted = false;

  /** The m bins. */
  protected IterMap<Chromosome, List<Integer>> mBins = new IterHashMap<Chromosome, List<Integer>>(25);

  /**
   * Instantiates a new binary gap search.
   */
  public BinaryGapSearch() {
    this(DEFAULT_BIN_SIZE);
  }

  /**
   * Instantiates a new binary gap search.
   *
   * @param binSize the bin size
   */
  public BinaryGapSearch(int binSize) {
    super(binSize);
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
    super.add(region, feature);

    mAutoSorted = false;
  }

  /**
   * Organize if not done so.
   */
  protected void organize() {
    if (!mAutoSorted) {

      mBins.clear();

      for (Entry<Chromosome, IterMap<Integer, GappedSearchFeatures<T>>> f : mFeatures) {
        mBins.put(f.getKey(), CollectionUtils.sortKeys(f.getValue()));
      }

      mAutoSorted = true;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#getFeatures(edu.
   * columbia.rdf.lib.bioinformatics.genome.Chromosome, int, int)
   */
  @Override
  public List<GappedSearchFeatures<T>> getFeatures(Chromosome chr, int start, int end) {
    // Make sure everything is sorted before doing anything
    organize();

    Map<Integer, GappedSearchFeatures<T>> features = mFeatures.get(chr);
    List<Integer> bins = mBins.get(chr);

    // SysUtils.err().println("bins", chr, mBins.keySet());

    if (features == null || features.size() == 0) {
      return Collections.emptyList();
    }

    int bs = start / mBinSize;
    int be = end / mBinSize;

    int is = getStartIndex(bins, bs);
    int ie = getEndIndex(bins, be);

    List<GappedSearchFeatures<T>> range = new ArrayList<GappedSearchFeatures<T>>();

    // System.err.println("is " + is + " " + ie + " " + start + " " + end + " "
    // + bs
    // + " " + be);
    // System.err.println("b " + bins.get(is) + " " + bins.get(ie));

    for (int i = is; i <= ie; ++i) {
      int bin = bins.get(i);

      GappedSearchFeatures<T> gsf = features.get(bin); // bins.get(i));

      if (gsf != null) {
        range.add(gsf);
      }
    }

    return range;
  }

  /**
   * public List<T> getClosestFeatures(GenomicRegion region, int n) { // Make sure
   * everything is sorted before doing anything organize();
   * 
   * Map<Integer, GappedSearchFeatures<T>> features =
   * mFeatures.get(region.getChr());
   * 
   * if (features == null) { return Collections.emptyList(); }
   * 
   * List<Integer> bins = mBins.get(region.getChr());
   * 
   * 
   * // The closest indices int is = getStartIndex(bins, region.getStart() /
   * mBinSize); int ie = getEndIndex(bins, region.getEnd() / mBinSize);
   * 
   * // find the index of the features closest to our point
   * 
   * int minD = Integer.MAX_VALUE; int closestIndex = Integer.MAX_VALUE;
   * 
   * for (int i = is; i <= ie; ++i) { int bin = bins.get(i);
   * 
   * int d = Math.abs(GenomicRegion.mid(region) -
   * features.get(bin).getPosition());
   * 
   * if (d < minD) { closestIndex = i; minD = d; } }
   * 
   * int closestP = features.get(bins.get(closestIndex)).getPosition();
   * 
   * // Sort the closest n points around this index to find 1st, 2nd, 3rd // etc
   * closest. Map<Integer, Integer> dMap = new HashMap<Integer, Integer>();
   * 
   * for (int i = 0; i <= n; ++i) { //System.err.println("c " + closestIndex + " "
   * + i + " " + n);
   * 
   * int it = closestIndex - i;
   * 
   * if (it >= 0) { dMap.put(Math.abs(features.get(bins.get(it)).getPosition() -
   * closestP), it); }
   * 
   * it = closestIndex + i;
   * 
   * if (it < features.size()) {
   * dMap.put(Math.abs(features.get(bins.get(it)).getPosition() - closestP), it);
   * } }
   * 
   * // If n = 0, thats the closest, n = 1, is the second closest int
   * closestIndexN = dMap.get(CollectionUtils.sort(dMap.keySet()).get(n));
   * 
   * GappedSearchFeatures<T> closestFeaturesN =
   * features.get(bins.get(closestIndexN));
   * 
   * List<T> ret = new UniqueArrayList<T>();
   * 
   * for (T item : closestFeaturesN) { ret.add(item); }
   * 
   * return ret; }
   *
   * @param <TT>         the generic type
   * @param closestIndex the closest index
   * @param features     the features
   * @param bins         the bins
   * @param used         the used
   * @param ret          the ret
   * @return the int
   */

  /**
   * Add features to the list of the nth closest features.
   * 
   * @param closestIndex
   * @param features
   * @param bins
   * @param used
   * @param ret
   * @return
   */
  private static final <TT> int addFeatures(int closestIndex, Map<Integer, GappedSearchFeatures<TT>> features,
      List<Integer> bins, Set<TT> used, List<List<TT>> ret) {
    int closestBin = bins.get(closestIndex);

    GappedSearchFeatures<TT> closestFeatures = features.get(closestBin);

    List<TT> l = new UniqueArrayList<TT>();

    for (Entry<GenomicRegion, List<TT>> r : closestFeatures) {
      for (TT item : r.getValue()) {
        if (!used.contains(item)) {
          l.add(item);
          used.add(item);
        }
      }
    }

    ret.add(l);

    return closestBin;
  }

  /**
   * Gets the start index in a list of ordered bins.
   *
   * @param bins  the bins
   * @param start the start
   * @return The start index in the list of bins that is closest (or contains) to
   *         the start.
   */
  public static int getStartIndex(List<Integer> bins, int start) {
    if (bins.size() < 2) {
      return 0;
    }

    // System.err.println("bins " + bins.get(0) + " " + bins.get(bins.size() -
    // 1) +
    // " " + bins.size() + " " + start);

    int is = 0;

    if (start <= bins.get(is)) {
      return is;
    }

    int ie = bins.size() - 1;

    if (start >= bins.get(ie)) {
      return ie;
    }

    // Binary search. Narrow the search to where point is lying between
    // two indices and pick the lower of the two

    int im = -1;
    int pm = -1;

    while (ie - is > 1) {
      im = (ie + is) / 2;

      pm = bins.get(im);

      // System.err.println("start " + start + " " + is + " " + ie + " " + im +
      // " " +
      // bins.get(is) + " " + bins.get(ie) + " " + pm);

      if (pm > start) {
        ie = im;
      } else if (pm < start) {
        is = im;
      } else {
        // We happen to exactly match the bin start
        return im; // bins.get(im); //return im;
      }
    }

    // If we made it this far, we have narrowed down the search to
    // two position so return the end

    /*
     * if (im == ie) { return ie; } else { return is; //ie; }
     */

    return is;
  }

  /**
   * Gets the closest index of features that either overlap this position or the
   * index of features just outside this point.
   *
   * @param bins the features
   * @param end  the end
   * @return the end index
   */
  public static int getEndIndex(List<Integer> bins, int end) {
    if (bins.size() < 2) {
      return bins.size() - 1;
    }

    int is = 0;

    if (end <= bins.get(is)) {
      return is;
    }

    int ie = bins.size() - 1;

    if (end >= bins.get(ie)) {
      return ie;
    }

    int im = -1;
    int pm = -1;

    while (ie - is > 1) {
      im = (ie + is) / 2;

      pm = bins.get(im);

      if (pm > end) {
        ie = im;
      } else if (pm < end) {
        is = im;
      } else {
        return im; // bins.get(im); //return im;
      }
    }

    // If we made it this far, we have narrowed down the search to
    // two position so return the end

    /*
     * if (im == is) { return is; } else { return ie; //ie; }
     */

    return ie;
  }

  /**
   * Contains chr.
   *
   * @param chr the chr
   * @return true, if successful
   */
  public boolean containsChr(Chromosome chr) {
    return mFeatures.containsKey(chr);
  }

  /**
   * Returns the total number of features for the chromosome.
   *
   * @param chr the chr
   * @return the int
   */
  public int size(Chromosome chr) {
    int ret = 0;

    IterMap<Integer, GappedSearchFeatures<T>> features = mFeatures.get(chr);

    for (Entry<Integer, GappedSearchFeatures<T>> f : features) {
      ret += f.getValue().size();
    }

    return ret;
  }

  /**
   * Gets the features at.
   *
   * @param chr the chr
   * @param i   the i
   * @return the features at
   */
  public GappedSearchFeatures<T> getFeaturesAt(Chromosome chr, int i) {
    return mFeatures.get(chr).get(mBins.get(chr).get(i));
  }

  /**
   * Gets the bins.
   *
   * @param chr the chr
   * @return the bins
   */
  public List<Integer> getBins(Chromosome chr) {
    return mBins.get(chr);
  }

  /**
   * Return the nth closest features by bin.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @param n     the n
   * @return the closest features
   */
  public List<List<T>> getClosestFeatures(Chromosome chr, int start, int end, int n) {
    // Make sure everything is sorted before doing anything
    organize();

    Map<Integer, GappedSearchFeatures<T>> features = mFeatures.get(chr);

    List<Integer> bins = mBins.get(chr);

    int bs = start;
    int be = end;

    if (mBinSize > 1) {
      bs /= mBinSize;
      be /= mBinSize;
    }

    // The closest indices
    int is = getStartIndex(bins, bs);
    int ie = getEndIndex(bins, be);

    // find the index of the features closest to our point

    int minD = Integer.MAX_VALUE;
    int closestIndex = Integer.MAX_VALUE;

    // System.err.println("is " + is + " " + ie);
    // System.err.println("huh " + features.keySet());

    for (int i = is; i <= ie; ++i) {
      int bin = bins.get(i);

      int d = Math.abs(GenomicRegion.mid(start, end) - bin);

      if (d < minD) {
        closestIndex = i;
        minD = d;
      }
    }

    List<List<T>> ret = new ArrayList<List<T>>(n);
    Set<T> used = new HashSet<T>();

    int closestBin = addFeatures(closestIndex, features, bins, used, ret);

    // SysUtils.err().println(chr, start, end, closestIndex);

    // int i = 0;
    int c = 0;

    int s = bins.size() - 1;

    int i1 = Math.min(s, closestIndex + 1);
    int i2 = Math.max(0, closestIndex - 1);

    int nthClosestIndex = closestIndex;

    while (c <= n) {
      int b1 = bins.get(i1);
      int b2 = bins.get(i2);

      // First we keep expanding the indices until we find genes
      // we haven't used before
      while (i1 < s) {
        boolean add = false;

        for (Entry<GenomicRegion, List<T>> r : features.get(b1)) {
          for (T item : r.getValue()) {
            if (!used.contains(item)) {
              add = true;
              break;
            }
          }
        }

        if (add) {
          break;
        }

        i1 = Math.min(s, i1 + 1);
        b1 = bins.get(i1);
      }

      while (i2 > 0) {
        boolean add = false;

        for (Entry<GenomicRegion, List<T>> r : features.get(b2)) {
          for (T item : r.getValue()) {
            if (!used.contains(item)) {
              add = true;
              break;
            }
          }
        }

        if (add) {
          break;
        }

        i2 = Math.max(0, i2 - 1);
        b2 = bins.get(i2);
      }

      // Once we have some bins to check, pick the closest on each
      // iteration.

      int d1 = Math.abs(b1 - closestBin);
      int d2 = Math.abs(b2 - closestBin);

      if (d1 <= d2) {
        nthClosestIndex = i1;

        // If the bin to the right of the current is closest then on
        // the next iteration, we need to check the bin to the left
        // with the next bin to the right. On the next iteration, the
        // bin on the left is probably the next closest, but we need
        // to check with the next closest on the right. We repeat
        // either incrementing i1 or decrementing i2 to find the
        // closest features and order them by absolute distance from
        // the current point.
        i1 = Math.min(s, i1 + 1);
      } else {
        nthClosestIndex = i2;
        i2 = Math.max(0, i2 - 1);
      }

      addFeatures(nthClosestIndex, features, bins, used, ret);

      ++c;
    }

    return ret;

    /*
     * int nthClosestBin = bins.get(nthClosestIndex);
     * 
     * GappedSearchFeatures<T> nthClosestFeatures = features.get(nthClosestBin);
     * 
     * List<T> ret = new UniqueArrayList<T>();
     * 
     * for (T item : nthClosestFeatures) { ret.add(item); }
     * 
     * return ret;
     */

    //
    // Look at at the closest items with a smaller start than the current
    // start.
    //

    /*
     * while (c <= n) { int it = closestIndex - i;
     * 
     * if (it < 0) { break; }
     * 
     * boolean add = true;
     * 
     * int bin = bins.get(it);
     * 
     * for (T item : features.get(bin)) { if (used.contains(item)) { add = false;
     * break; } }
     * 
     * // Essentially look for the closest unique genes ranked by distance // We do
     * not want the same gene listed x times because it // occupies multiple bins if
     * (add) { // Log the absolute position
     * dMap.put(Math.abs(features.get(bin).getPosition() - closestP), it);
     * 
     * for (T item : features.get(bin)) { used.add(item); }
     * 
     * ++c; }
     * 
     * ++i; }
     * 
     * // // Repeat looking at positions greater than the start //
     * 
     * i = 0; c = 0;
     * 
     * while (c <= n) { //System.err.println("c " + closestIndex + " " + i + " " +
     * n);
     * 
     * int it = closestIndex + i;
     * 
     * if (it >= features.size()) { break; }
     * 
     * boolean add = true;
     * 
     * int bin = bins.get(it);
     * 
     * for (T item : features.get(bin)) { if (used.contains(item)) { add = false;
     * break; } }
     * 
     * if (add) { dMap.put(Math.abs(features.get(bin).getPosition() - closestP),
     * it);
     * 
     * for (T item : features.get(bin)) { used.add(item); }
     * 
     * ++c; }
     * 
     * ++i; }
     * 
     * // Sort the distances from the position so we now have an ordering // of
     * positions from 0 - (n - 1)th closest List<Integer> sortedDistances =
     * CollectionUtils.sortKeys(dMap);
     * 
     * // If n = 0, thats the closest, n = 1, is the second closest etc int
     * closestIndexN = sortedDistances.get(Mathematics.bound(n, 0, dMap.size() -
     * 1));
     * 
     * int bin = bins.get(closestIndexN);
     * 
     * GappedSearchFeatures<T> closestFeaturesN = features.get(bin);
     * 
     * List<T> ret = new UniqueArrayList<T>();
     * 
     * for (T item : closestFeaturesN) { ret.add(item); }
     * 
     * return ret;
     */
  }

  /**
   * Get the closest features distance n from location. For example if n = 0,
   * return the closest, n = 1, return the second closest, n = 3 the third closest
   * etc.
   *
   * @param region the region
   * @param n      the n
   * @return the closest features
   */
  public List<List<T>> getClosestFeatures(GenomicRegion region, int n) {
    return getClosestFeatures(region.getChr(), region.getStart(), region.getEnd(), n);
  }

}
