package org.jebtk.bioinformatics.genomic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map.Entry;

import org.jebtk.core.search.SearchQuery;

public class GeneSearchQuery extends SearchQuery<GenomicElement> {

  private GFBGenes mGenes;
  private GenomicType mType;

  public GeneSearchQuery(GFBGenes genes) {
    this(genes, GenomicType.GENE);
  }

  public GeneSearchQuery(GFBGenes genes, GenomicType type) {
    mGenes = genes;
    mType = type;
  }

  @Override
  public Collection<GenomicElement> match(String s, boolean exact, boolean include) {

    String ls = s.toLowerCase();

    List<GenomicElement> ret = new ArrayList<GenomicElement>();

    // System.err.println("qs " + s);

    List<GenomicElement> genes = mGenes.getGenes(ls, mType);

    // if (!include) {
    // Add all the genes not in this list of genes
    // genes = CollectionUtils.compliment(mGenes.getGenes(), genes);
    // }

    if (exact) {
      genes = exact(ls, genes);
    }

    ret.addAll(genes);

    return ret;
  }

  /**
   * Return only the exact matches.
   * 
   * @param s
   * @param genes
   * @return
   */
  private static List<GenomicElement> exact(String s, Collection<GenomicElement> genes) {
    List<GenomicElement> ret = new ArrayList<GenomicElement>(genes.size());

    for (GenomicElement gene : genes) {
      boolean found = false;

      for (Entry<String, Object> e : gene.getProperties()) {
        if (e.getValue().toString().toLowerCase().equals(s)) {
          found = true;
          break;
        }
      }

      if (found) {
        ret.add(gene);
      }
    }

    return ret;
  }

}
