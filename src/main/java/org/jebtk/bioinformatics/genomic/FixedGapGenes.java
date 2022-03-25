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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jebtk.bioinformatics.gapsearch.FixedGapSearch;
import org.jebtk.bioinformatics.gapsearch.GappedSearchFeatures;
import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.io.PathUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public class FixedGapGenes extends SingleGenesDB {
  // private static final GeneService INSTANCE = new GeneService();

  private static final long serialVersionUID = 1L;

  /**
   * The constant DEFAULT_RES.
   */
  public static final String DEFAULT_RES = "res/ucsc_refseq_exons_entrez_hg19.txt.gz";

  /**
   * The constant DEFAULT_FILE.
   */
  public static final Path DEFAULT_FILE = PathUtils.getPath(DEFAULT_RES);

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(FixedGapGenes.class);

  private FixedGapSearch<GenomicElement> mSearch = new FixedGapSearch<GenomicElement>();

  /**
   * The member symbol main variant.
   */
  // private Map<String, GenomicElement> mSymbolMainVariant = new HashMap<String,
  // GenomicElement>();

  /**
   * The member entrez map.
   */
  // private Map<String, Set<Gene>> mEntrezMap =
  // DefaultHashMap.create(new HashSetCreator<Gene>());

  /**
   * The member symbol map.
   */

  // private IterMap<String, Set<String>> mTypeMap = DefaultHashMap
  // .create(new TreeSetCreator<String>());

  // private IterMap<String, List<GenomicElement>> mIdMap = DefaultHashMap
  // .create(new UniqueArrayListCreator<GenomicElement>());

  /**
   * The member ref seq map.
   */
  // private Map<String, Gene> mRefSeqMap =
  // new HashMap<String, Gene>();

  public FixedGapGenes(Genome genome) {
    super(genome);
  }

  @Override
  public void add(GenomicElement gene) {
    mSearch.add(gene, gene);

    // Map to exons

    /*
     * if (gene.getExonCount() > 0) { for (GenomicRegion exon : gene) {
     * super.add(exon, gene); } } else { super.add(region, gene); }
     */

    // for (String name : gene.getPropertyNames()) {
    // String value = gene.getProperty(name);
    // mTypeMap.get(name).add(value);
    // mIdMap.get(sanitize(value)).add(gene);
    // }
  }

//  /**
//   * Auto find main variants.
//   */
//  public void autoFindMainVariants() {
//
//    if (mFindMainVariants) {
//      //
//      // Find the representative gene e.g. variant 1
//      //
//
//      Set<String> symbols = mTypeMap.get(Gene.GENE_NAME_TYPE);
//
//      for (String name : symbols) {
//        List<GenomicElement> genes = mIdMap.get(sanitize(name));
//
//        // try and find variant 1
//
//        int maxWidth = 0;
//
//        for (GenomicElement gene : genes) {
//          int width = gene.mLength;
//
//          if (width > maxWidth) {
//            mSymbolMainVariant.put(name, gene);
//
//            maxWidth = width;
//          }
//        }
//      }
//
//      mFindMainVariants = false;
//    }
//  }

  /*
   * public GenomicElement lookup(String genome, String id) throws IOException {
   * GenomicElement gene = getGene(genome, id);
   * 
   * if (gene != null) { return gene; }
   * 
   * gene = findMainVariant(id);
   * 
   * return gene; }
   */

//  @Override
//  public List<GenomicElement> getElements(Genome genome,
//      String id,
//      GenomicType type) throws IOException {
//    return mIdMap.get(sanitize(id));
//  }

  @Override
  public List<GenomicElement> getElements() {
    return mSearch.getFeatures();
  }

  /**
   * Find genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException
   */
  @Override
  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
    List<GenomicElement> features = mSearch.find(region, minBp);

    return filterByType(features, type);
  }

  /**
   * Find closest genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException
   */
  @Override
  public List<GenomicElement> closest(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    List<GenomicElement> features = mSearch.getClosestFeatures(region);

    return filterByType(features, type);
  }

  @Override
  public List<List<GenomicElement>> nthClosest(Genome genome, GenomicRegion region, int n, GenomicType type) {
    // return mSearch.getClosestFeatures(region, n);

    return nthClosest(region.getChr(), region.getStart(), region.getEnd(), n, type);
  }

  private List<List<GenomicElement>> nthClosest(Chromosome chr, int start, int end, int n, GenomicType type) {
    int bs = mSearch.getBin(start);
    int be = mSearch.getBin(end);

    Map<Integer, GappedSearchFeatures<GenomicElement>> features = mSearch.get(chr);

    if (features.size() == 0) {
      return Collections.emptyList();
    }

    List<Integer> bins = CollectionUtils.sortKeys(features);

    int is = bins.indexOf(bs);
    int ie = bins.indexOf(be);

    int bn = bins.size();
    int l = bins.size() - 1;

    List<GappedSearchFeatures<GenomicElement>> bf = null;

    // SysUtils.err().println(bs, be, is, ie);

    int mid = (start + end) / 2;

    do {
      is = Mathematics.bound(is, 0, l);
      ie = Mathematics.bound(ie, 0, l);

      bs = bins.get(is);
      be = bins.get(ie);

      // Keep expanding bin search area around location until we find enough
      // items to order by 1st, 2nd, 3rd... closest.

      bf = mSearch.getFeaturesByBin(chr, bs, be);

      //
      // Count
      //

      IterMap<Integer, List<GenomicElement>> closestMap = DefaultTreeMap.create(new ArrayListCreator<GenomicElement>());

      for (GappedSearchFeatures<GenomicElement> gsf : bf) {
        for (Entry<GenomicRegion, List<GenomicElement>> r : gsf) {
          // distance from item to

          int d;
          GenomicRegion region = r.getKey();

          if (Strand.isSense(region.getStrand())) {
            d = Math.abs(region.getStart() - mid);
          } else {
            d = Math.abs(region.getEnd() - mid);
          }

          for (GenomicElement item : r.getValue()) {
            if (item.mType.equals(type)) {
              closestMap.get(d).add(item);
            }

            // SysUtils.err().println("closest", n, d, item, mid,
            // region.getStart());
          }
        }
      }

      if (closestMap.size() >= n) {
        // System.err.println("closest" + closestMap);
        // Once we have enough closest genes, assemble and return
        List<List<GenomicElement>> ret = new ArrayList<List<GenomicElement>>(n);

        int c = 0;

        for (Entry<Integer, List<GenomicElement>> entry : closestMap) {
          ret.add(entry.getValue());

          ++c;

          if (c == n) {
            break;
          }
        }

        return ret;
      }

      //
      //
      //

      // Expand search region
      --is;
      ++ie;
    } while (is > -1 || ie < bn);

    return Collections.emptyList();
  }

  /*
   * @Override public GenomicElement findMainVariant(String name) {
   * autoFindMainVariants();
   * 
   * return mSymbolMainVariant.get(name.toUpperCase()); }
   */

//  @Override
//  public Iterable<String> getIds(String type) {
//    return mTypeMap.get(type);
//  }

  @Override
  public boolean contains(Chromosome chr) {
    return mSearch.contains(chr);
  }

}