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
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.jebtk.core.event.ChangeListeners;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public abstract class GenomicElementsDB extends ChangeListeners {

  private static final long serialVersionUID = 1L;

  /** Empty gene set that can be used as a placeholder */
  public static final GenomicElementsDB EMPTY = new GenomicElementsDB() {

    private static final long serialVersionUID = 1L;

    @Override
    public Iterable<Genome> getGenomes() {
      return Collections.emptyList();
    }

    @Override
    public List<GenomicElement> getElements(Genome g, String search, GenomicType type) {
      return Collections.emptyList();
    }

    @Override
    public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
      return Collections.emptyList();
    }

    @Override
    public void add(GenomicElement element) {
      // Do nothing
    }
  };

  public abstract void add(GenomicElement element);

  /**
   * Lookup a gene by either symbol or refseq.
   *
   * @param id the id
   * @return the gene
   * @throws IOException
   */
  // public abstract GenomicElement lookup(String genome, String id) throws
  // IOException;

  public GenomicElement getElement(Genome genome, String search, GenomicType type) throws IOException {
    System.err.println("Searching for genes in " + genome + " " + this.getClass());

    List<GenomicElement> genes = getElements(genome, search, type);

    if (genes.size() > 0) {
      return genes.get(0);
    } else {
      return null;
    }
  }

  public List<GenomicElement> getElements(Genome g, String search, GenomicType type) {
    return Collections.emptyList();
  }

  public List<GenomicElement> getElements() {
    return Collections.emptyList();
  }

  public List<GenomicElement> getElements(Chromosome chr) {
    return Collections.emptyList();
  }

  public List<GenomicElement> find(GenomicRegion region) {
    return find(Genome.NO_GENOME, region);
  }

  public List<GenomicElement> find(Genome genome, GenomicRegion region) {
    return find(genome, region, GenomicType.REGION);
  }

  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type) {
    return find(genome, region, type, 1);
  }

  /**
   * Find genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException
   */
  public abstract List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp);

  /**
   * Find closest genes.
   *
   * @param region the region
   * @return the list
   * @throws IOException
   */
  public List<GenomicElement> closest(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {

    List<GenomicElement> elements = find(genome, region, type, minBp);

    return closest(region, elements);
  }

  public List<List<GenomicElement>> nthClosest(Genome genome, GenomicRegion region, int n, GenomicType type) {
    return Collections.emptyList();
  }

  /**
   * Find closest genes by tss.
   *
   * @param region the region
   * @return the list
   * @throws IOException
   */
  public List<GenomicElement> closestByTss(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    Collection<GenomicElement> genes = closest(genome, region, type, minBp); // findGenes(region);

    List<GenomicElement> ret = new ArrayList<GenomicElement>();

    int minD = Integer.MAX_VALUE;

    for (GenomicElement gene : genes) {
      GenomicRegion tss = Gene.tssRegion(gene);

      minD = Math.min(minD, GenomicRegion.midAbsDist(region, tss));
    }

    for (GenomicElement gene : genes) {
      GenomicRegion tss = Gene.tssRegion(gene);

      int d = GenomicRegion.midAbsDist(region, tss);

      if (d == minD) {
        ret.add(gene);
      }
    }

    return ret;
  }

  /*
   * public List<GenomicElement> search(GenomicRegion region) {
   * List<GenomicElement> elements = GenomicRegions
   * .getFixedGapSearch(getElements(region.mChr)) .getValues(region);
   * 
   * return elements; }
   */

  /**
   * Return the set of ids (e.g. RefSeq ids) associated with a given id type.
   *
   * @param type the type
   * @return the ids
   */
  public Iterable<String> getIds(String type) {
    return Collections.emptyList();
  }

  public boolean contains(Chromosome chr) {
    return false;
  }

  /**
   * Should return the supported databases
   * 
   * @return
   */
  public Iterable<Genome> getGenomes() {
    return Collections.emptyList();
  }

  //
  // Static methods
  //

  public static String sanitize(String name) {
    return name.toUpperCase();
  }

  public static GFF3Parser gff3Parser() {
    return new GFF3Parser();
  }

  public static GTB1Parser gtbParser() {
    return new GTB1Parser();
  }

  public static GTB2Parser gtb2Parser() {
    return new GTB2Parser();
  }

  public List<GenomicElement> overlapping(Genome genome, GenomicRegion region, GenomicType type) throws IOException {
    return overlapping(genome, region, type, 1);
  }

  public List<GenomicElement> overlapping(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    List<GenomicElement> ret = new ArrayList<GenomicElement>();

    overlapping(genome, region, type, minBp, ret);

    return ret;
  }

  public void overlapping(Genome genome, GenomicRegion region, GenomicType type, int minBp, List<GenomicElement> ret)
      throws IOException {
    List<GenomicElement> elements = find(genome, region, type, minBp);

    if (elements.size() == 0) {
      return;
    }

    for (GenomicElement g : elements) {
      GenomicRegion overlap = GenomicRegion.overlap(region, g);

      // SysUtils.err().println("overlap", overlap.getLength(), region, g);

      if (overlap != null && overlap.getLength() >= minBp) {
        ret.add(g);
      }
    }
  }

  private static List<GenomicElement> closest(GenomicRegion region, List<GenomicElement> elements) {

    List<GenomicElement> ret = new ArrayList<GenomicElement>(elements.size());

    int mid = GenomicRegion.mid(region);
    int minD = Integer.MAX_VALUE;

    for (GenomicElement element : elements) {
      int d = Math.abs(mid - element.getStart()); // GenomicRegion.mid(gene));

      if (d < minD) {
        minD = d;
      }
    }

    for (GenomicElement element : elements) {
      int d = Math.abs(mid - element.getStart()); // GenomicRegion.mid(gene));

      if (d == minD) {
        ret.add(element);
      }
    }

    return ret;
  }
}