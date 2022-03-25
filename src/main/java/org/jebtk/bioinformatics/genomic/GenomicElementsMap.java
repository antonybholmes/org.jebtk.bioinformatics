package org.jebtk.bioinformatics.genomic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jebtk.bioinformatics.gapsearch.GapSearch;
import org.jebtk.core.SizeGetter;
import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;

public class GenomicElementsMap extends GenomicElementsDB
    implements SizeGetter, Iterable<Entry<Chromosome, List<GenomicElement>>> {

  private static final long serialVersionUID = 1L;

  private IterMap<Chromosome, List<GenomicElement>> mElementMap = DefaultTreeMap
      .create(new ArrayListCreator<GenomicElement>());

  private Map<Chromosome, GapSearch<GenomicElement>> mFindMap = new HashMap<Chromosome, GapSearch<GenomicElement>>();

  private int mSize = 0;

  @Override
  public Iterator<Entry<Chromosome, List<GenomicElement>>> iterator() {
    return mElementMap.iterator();
  }

  @Override
  public List<GenomicElement> getElements(Chromosome chr) {
    return Collections.unmodifiableList(mElementMap.get(chr));
  }

  @Override
  public List<GenomicElement> getElements() {
    List<GenomicElement> ret = new ArrayList<GenomicElement>();

    for (Entry<Chromosome, List<GenomicElement>> item : this) {
      ret.addAll(item.getValue());
    }

    return ret;
  }

  public int size() {
    return mSize;
  }

  @Override
  public void add(GenomicElement e) {
    mElementMap.get(e.mChr).add(e);

    mFindMap.remove(e.mChr);

    ++mSize;
  }

  @Override
  public List<GenomicElement> getElements(Genome g, String search, GenomicType type) {
    return Collections.emptyList();
  }

  @Override
  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
    if (!mFindMap.containsKey(region.mChr)) {
      mFindMap.put(region.mChr, GenomicRegions.getFixedGapSearch(mElementMap.get(region.mChr)));
    }
    // return GenomicRegions.getFixedGapSearch(mElementMap.get(region.mChr))
    // .getValues(region);

    return mFindMap.get(region.mChr).find(region, minBp);
  }

  public void addAll(Collection<GenomicElement> elements) {
    for (GenomicElement e : elements) {
      add(e);
    }

  }

}
