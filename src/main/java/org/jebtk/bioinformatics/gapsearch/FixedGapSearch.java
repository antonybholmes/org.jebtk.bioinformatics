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
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.TreeMapCreator;
import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.sys.SysUtils;

/**
 * Use fixed size blocks to find features.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class FixedGapSearch<T> extends GapSearch<T> {

  /**
   * The constant DEFAULT_BIN_SIZE.
   */
  protected static final int DEFAULT_BIN_SIZE = 10000;

  /**
   * The member features.
   */
  protected IterMap<Chromosome, IterMap<Integer, GappedSearchFeatures<T>>> mFeatures = DefaultTreeMap
      .create(new TreeMapCreator<Integer, GappedSearchFeatures<T>>());

  // protected IterMap<Chromosome, BinMap<GappedSearchFeatures<T>>> mFeatures;

  /**
   * The member size.
   */
  protected int mSize = 0;

  protected final int mBinSize;

  /**
   * The member bin size.
   */
  // protected final int mBinSize;

  /**
   * Instantiates a new fixed gap search.
   */
  public FixedGapSearch() {
    this(DEFAULT_BIN_SIZE);
  }

  /**
   * Instantiates a new fixed gap search.
   *
   * @param binSize the bin size
   */
  public FixedGapSearch(int binSize) {
    mBinSize = Math.max(1, binSize);

    /*
     * mFeatures = DefaultTreeMap.create(new
     * EntryCreator<BinMap<GappedSearchFeatures<T>>>(){
     * 
     * @Override public BinMap<GappedSearchFeatures<T>> newEntry() { return new
     * BinMap<GappedSearchFeatures<T>>(); }});
     */
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
    int startBin = getBin(region.getStart());
    int endBin = getBin(region.getEnd());

    Chromosome chr = region.getChr();

    for (int bin = startBin; bin <= endBin; ++bin) {
      if (!mFeatures.get(chr).containsKey(bin)) {
        mFeatures.get(chr).put(bin, new GappedSearchFeatures<T>(bin));
      }

      mFeatures.get(chr).get(bin).add(region, feature);
    }

    ++mSize;
  }

  public void addAll(GenomicRegion region, Collection<T> features) {
    int startBin = getBin(region.getStart());
    int endBin = getBin(region.getEnd());

    Chromosome chr = region.getChr();

    for (int bin = startBin; bin <= endBin; ++bin) {
      if (!mFeatures.get(chr).containsKey(bin)) {
        mFeatures.get(chr).put(bin, new GappedSearchFeatures<T>(bin));
      }

      for (T feature : features) {
        mFeatures.get(chr).get(bin).add(region, feature);
      }
    }

    ++mSize;
  }

  public int getBin(int x) {
    return x / mBinSize;
  }

  /**
   * Adds the feature.
   *
   * @param chr     the chr
   * @param bin     the bin
   * @param feature the feature
   */
  // public void addFeature(Chromosome chr, int bin, T feature) {
  // if (!mFeatures.get(chr).containsKey(bin)) {
  // mFeatures.get(chr).put(bin, new GappedSearchFeatures<T>(bin));
  // }
  //
  // mFeatures.get(chr).get(bin).add(feature);
  //
  // ++mSize;
  // }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#size()
   */
  @Override
  public int size() {
    return mSize;
  }

  @Override
  public boolean contains(Chromosome chr) {
    return mFeatures.containsKey(chr);
  }

  public IterMap<Integer, GappedSearchFeatures<T>> get(Chromosome chr) {
    return mFeatures.get(chr);
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#getFeatureList()
   */
  @Override
  public List<T> getFeatures() {
    List<T> ret = new UniqueArrayList<T>();

    for (Entry<Chromosome, IterMap<Integer, GappedSearchFeatures<T>>> f : mFeatures) {
      IterMap<Integer, GappedSearchFeatures<T>> starts = f.getValue();

      for (Entry<Integer, GappedSearchFeatures<T>> e : starts) {
        GappedSearchFeatures<T> sf = e.getValue();

        for (Entry<GenomicRegion, List<T>> r : sf) {
          ret.addAll(r.getValue());
        }
      }

      /*
       * for (int start : mFeatures.get(chr)) { for (GenomicRegion region :
       * mFeatures.get(chr).get(start)) { for (T item :
       * mFeatures.get(chr).get(start).getValues(region)) { ret.add(item); } } }
       */
    }

    return ret;
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#getFeatureList(edu.
   * columbia.rdf.lib.bioinformatics.genome.Chromosome)
   */
  @Override
  public List<T> getFeatures(Chromosome chr) {
    List<T> ret = new UniqueArrayList<T>();

    for (Entry<Integer, GappedSearchFeatures<T>> f : mFeatures.get(chr)) {
      for (Entry<GenomicRegion, List<T>> r : f.getValue()) {
        ret.addAll(r.getValue());
      }
    }

    return ret;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.gapsearch.GapSearch#getFeatures(edu.
   * columbia.rdf.lib.bioinformatics.genome.Chromosome, int, int)
   */
  @Override
  public List<GappedSearchFeatures<T>> getFeatures(Chromosome chr, int start, int end) {
    int is = getBin(start);
    int ie = getBin(end);

    return getFeaturesByBin(chr, is, ie);
  }

  /**
   * Gets the features by bin.
   *
   * @param chr  the chr
   * @param sbin the sbin
   * @param ebin the ebin
   * @return the features by bin
   */
  public List<GappedSearchFeatures<T>> getFeaturesByBin(Chromosome chr, int sbin, int ebin) {
    Map<Integer, GappedSearchFeatures<T>> features = mFeatures.get(chr);

    if (features.size() == 0) {
      return Collections.emptyList();
    }

    List<GappedSearchFeatures<T>> range = new ArrayList<GappedSearchFeatures<T>>();

    for (int i = sbin; i <= ebin; ++i) {
      GappedSearchFeatures<T> f = features.get(i);

      if (f != null) {
        range.add(f);
      }
    }

    return range;
  }

  /**
   * Return the nth closest features.
   * 
   * @param region
   * @param n
   * @return
   */
  public List<List<T>> getClosestFeatures(GenomicRegion region, int n) {
    return getClosestFeatures(region.getChr(), region.getStart(), region.getEnd(), n);
  }

  public List<List<T>> getClosestFeatures(Chromosome chr, int start, int end, int n) {
    int bs = getBin(start);
    int be = getBin(end);

    Map<Integer, GappedSearchFeatures<T>> features = mFeatures.get(chr);

    List<Integer> bins = CollectionUtils.sortKeys(features);

    int is = bins.indexOf(bs);
    int ie = bins.indexOf(be);

    int l = bins.size() - 1;

    List<GappedSearchFeatures<T>> bf = null;

    SysUtils.err().println(bs, be, is, ie);

    int mid = (start + end) / 2;

    while (is >= 0 || ie < bins.size()) {
      bs = bins.get(is);
      be = bins.get(ie);

      // Keep expanding bin search area around location until we find enough
      // items to order by 1st, 2nd, 3rd... closest.

      bf = getFeaturesByBin(chr, bs, be);

      //
      // Count
      //

      IterMap<Integer, List<T>> closestMap = DefaultTreeMap.create(new ArrayListCreator<T>());

      for (GappedSearchFeatures<T> gsf : bf) {
        for (Entry<GenomicRegion, List<T>> r : gsf) {
          // distance from item to
          int d;
          GenomicRegion region = r.getKey();

          if (Strand.isSense(region.getStrand())) {
            d = Math.abs(region.getStart() - mid);
          } else {
            d = Math.abs(region.getEnd() - mid);
          }

          closestMap.get(d).addAll(r.getValue());

          if (closestMap.size() == n) {
            // Once we have enough closest genes, assemble and return
            List<List<T>> ret = new ArrayList<List<T>>(n);

            for (Entry<Integer, List<T>> entry : closestMap) {
              ret.add(entry.getValue());
            }

            return ret;
          }
        }
      }

      //
      //
      //

      // Expand search region
      is = Math.max(0, is - 1);
      ie = Math.min(ie + 1, l);
    }

    return Collections.emptyList();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<Chromosome> iterator() {
    return mFeatures.keySet().iterator();
  }

  /**
   * Adds a search object to this one.
   *
   * @param gappedSearch the gapped search
   */
  public void add(FixedGapSearch<T> gappedSearch) {
    for (Entry<Chromosome, IterMap<Integer, GappedSearchFeatures<T>>> f1 : gappedSearch.mFeatures) {
      for (Entry<Integer, GappedSearchFeatures<T>> f2 : f1.getValue()) {
        GappedSearchFeatures<T> search = f2.getValue();

        for (Entry<GenomicRegion, List<T>> r : search) {
          addAll(r.getKey(), r.getValue());
        }
      }

      /*
       * for (int bin : gappedSearch.mFeatures.get(chr)) {
       * 
       * 
       * for (GenomicRegion region : gappedSearch.mFeatures.get(chr).get(bin)) { for
       * (T item : mFeatures.get(chr).get(bin).getValues(region)) { add(region, item);
       * } } }
       */
    }
  }

  /**
   * Return the closest set of features to the mid-point of a region.
   *
   * @param region the region
   * @return the closest features
   */
  @Override
  public List<T> getClosestFeatures(GenomicRegion region) {
    if (region == null) {
      return Collections.emptyList();
    }

    // organize();

    List<GappedSearchFeatures<T>> allFeatures = getFeatures(region);

    if (allFeatures.size() == 0) {
      return Collections.emptyList();
    }

    int minD = Integer.MAX_VALUE;

    GappedSearchFeatures<T> minF = null;

    int mid = GenomicRegion.mid(region);

    for (GappedSearchFeatures<T> features : allFeatures) {
      int d = Math.abs(mid - features.getPosition());

      if (d < minD) {
        minF = features;
        minD = d;
      }
    }

    List<T> ret = new UniqueArrayList<T>();

    for (Entry<GenomicRegion, List<T>> r : minF) {
      ret.addAll(r.getValue());
    }

    return ret;
  }

}
