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

import java.awt.Color;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.jebtk.core.collections.DefaultHashMap;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterHashMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.TreeSetCreator;
import org.jebtk.core.collections.UniqueArrayListCreator;
import org.jebtk.core.text.TextUtils;

import com.fasterxml.jackson.annotation.JsonGetter;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;

// TODO: Auto-generated Javadoc
/**
 * Represents a genomic entity such as a gene or exon. An entity is a genomic
 * region with a collection of annotations and sub entities to describe a
 * genomic feature such as a gene, transcript or exon. It is designed for use
 * with different annotation databases such as Refseq and Ensembl where
 * different ids and nomenclature are used.
 */
@JsonPropertyOrder({ "loc", "strand", "type", "properties", "tags" })
public class GenomicElement extends GenomicRegion {
  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /** The m id map. */
  protected IterMap<String, Object> mPropertyMap = new IterHashMap<String, Object>();

  /** The m tags. */
  private Set<String> mTags = new TreeSet<String>();

  /** The m elem map. */
  private IterMap<GenomicType, List<GenomicElement>> mElemMap = DefaultHashMap
      .create(new UniqueArrayListCreator<GenomicElement>());

  /** The m utr 5 p. */
  // private List<Exon> mUtr5p = new ArrayList<Exon>();

  /** The m text. */
  private String mText;

  /** The m parent. */
  private GenomicElement mParent = null;

  /** The m type. */
  @JsonIgnore
  public final GenomicType mType;

  private Color mColor = Color.BLACK;

  /**
   * Instantiates a new genomic entity.
   *
   * @param type the type
   * @param l    the l
   */
  public GenomicElement(GenomicType type, GenomicRegion l) {
    super(l);

    mType = type;
  }

  /**
   * Instantiates a new genomic entity.
   *
   * @param type the type
   * @param l    the l
   * @param s    the s
   */
  public GenomicElement(GenomicType type, GenomicRegion l, Strand s) {
    super(l, s);

    mType = type;
  }

  public GenomicElement(GenomicType type, Chromosome chr, int start, int end) {
    this(type, chr, start, end, Strand.SENSE);
  }

  public GenomicElement(GenomicType type, Chromosome chr, int start, int end, Strand strand) {
    super(chr, start, end, strand);

    mType = type;
  }

  /**
   * Gets the type.
   *
   * @return the type
   */
  @JsonGetter("type")
  public GenomicType getType() {
    return mType;
  }

  @JsonIgnore
  public Color getColor() {
    return mColor;
  }

  public GenomicElement setColor(Color color) {
    return setColor(color, true);
  }

  public GenomicElement setColor(Color color, boolean propogate) {
    setColor(this, color, propogate);

    return this;
  }

  /**
   * Set a parent entity to indicate this is a child of another.
   *
   * @param gene the new parent
   */
  public void setParent(GenomicElement gene) {
    mParent = gene;
  }

  /**
   * Gets the parent.
   *
   * @return the parent
   */
  @JsonIgnore
  public GenomicElement getParent() {
    return mParent;
  }

  /**
   * Add an entity as a child of this one. For example an exon could be a child of
   * transcript.
   * 
   * @param e A genomic entity.
   */
  public void addChild(GenomicElement e) {
    mElemMap.get(e.mType).add(e);
  }

  /**
   * Gets the child types.
   *
   * @return the child types
   */
  @JsonIgnore
  public Iterable<GenomicType> getChildTypes() {
    return mElemMap.keySet().stream().sorted().collect(Collectors.toList());
  }

  /**
   * Return child entities of this one of a specific type.
   *
   * @param type a genomic type.
   * @return the children
   */
  @JsonIgnore
  public Iterable<GenomicElement> getChildren(GenomicType type) {
    return mElemMap.get(type);
  }

  @JsonIgnore
  public Set<Entry<GenomicType, List<GenomicElement>>> getChildren() {
    return mElemMap.entrySet();
  }

  /**
   * Gets the child count.
   *
   * @param type the type
   * @return the child count
   */
  @JsonIgnore
  public int getChildCount(GenomicType type) {
    return mElemMap.get(type).size();
  }

  // public Iterable<Entry<String, String>> getPropertyNames() {
  // return mPropertyMap;
  // }

  /**
   * Sets the id.
   *
   * @param type the type
   * @param name the name
   * @return the genomic entity
   */
  public GenomicElement setProperty(String name, String value) {
    mPropertyMap.put(name, value);

    return this;
  }

  public GenomicElement setProperty(String name, int value) {
    mPropertyMap.put(name, value);

    return this;
  }

  public GenomicElement setProperty(String name, double value) {
    mPropertyMap.put(name, value);

    return this;
  }

  /*
   * public GenomicElement setProperty(String name, Tag property) { if
   * (TextUtils.isNullOrEmpty(name)) { return this; }
   * 
   * name = Tag.format(name);
   * 
   * mPropertyMap.put(name, property);
   * 
   * return this; }
   */

  /**
   * Returns true if the entity contains an id with a given name.
   * 
   * @param name the name of the id.
   * @return true if the name exists, false otherwise.
   */
  @JsonIgnore
  public boolean hasProperty(String name) {
    return mPropertyMap.containsKey(name);
  }

  /**
   * Return the different types of ids associated with this entity. Typically this
   * will be 'refseq' or 'gene_symbol'.
   *
   * @return the ids
   */
  @JsonIgnore
  public Iterable<String> getPropertyNames() {
    return mPropertyMap.entrySet().stream().sorted(Map.Entry.comparingByKey()).map(e -> e.getKey())
        .collect(Collectors.toList());
  }

  /**
   * Gets the id count.
   *
   * @return the id count
   */
  @JsonIgnore
  public int getPropertyCount() {
    return mPropertyMap.size();
  }

  /**
   * Gets the id.
   *
   * @param type the type
   * @return the id
   */
  @JsonIgnore
  public String getProperty(GeneIdType type) {
    return getProperty(type.toString());
  }

  /**
   * Return a property. If the property does not exist, a property with the string
   * value of 'n/a' will be automatically created so that a null is never
   * returned.
   *
   * @param type the type
   * @return the id
   */
  private Object _getProperty(String name) {
    name = Tag.format(name);

    // System.err.println("tag:" + name);

    if (!mPropertyMap.containsKey(name)) {
      setProperty(name, TextUtils.NA);
    }

    return mPropertyMap.get(name);
  }

  public String getProperty(String name) {
    return _getProperty(name).toString();
  }

  public int getInt(String name) {
    return (Integer) _getProperty(name);
  }

  public double getDouble(String name) {
    return (Double) _getProperty(name);
  }

  /**
   * Add an arbitrary tag to the entity such as an some meta data better
   * describing it.
   *
   * @param tag the tag
   * @return
   */
  public GenomicElement addTag(String tag) {
    mTags.add(tag);

    return this;
  }

  public GenomicElement addTag(int tag) {
    return addTag(Integer.toString(tag));
  }

  public GenomicElement addTag(double tag) {
    return addTag(Double.toString(tag));
  }

  public GenomicElement addTags(Collection<String> tags) {
    mTags.addAll(tags);

    return this;
  }

  /**
   * Gets the tags.
   *
   * @return the tags
   */
  @JsonGetter("tags")
  public Iterable<String> getTags() {
    return mTags;
  }

  /**
   * Gets the tag count.
   *
   * @return the tag count
   */
  @JsonIgnore
  public int getTagCount() {
    return mTags.size();
  }

  /**
   * Returns the element with the start and end based on orientation so that the
   * start may be after the end
   *
   * @return the tss
   */
  @JsonIgnore
  public GenomicRegion getTss() {
    return tssRegion(this);
  }

  /*
   * @Override public boolean equals(Object o) { if (o instanceof GenomicEntity) {
   * return toString().equals(((GenomicEntity) o).toString()); } else { return
   * super.equals(o); } }
   */

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genomic.GenomicRegion#compareTo(org.jebtk.
   * bioinformatics.genomic.Region)
   */
  @Override
  public int compareTo(Region r) {
    int c = super.compareTo(r);

    // Different location so sort
    if (c != 0) {
      return c;
    }

    // Same location, have to test toString

    if (r instanceof GenomicElement) {
      c = toString().compareTo(((GenomicElement) r).toString());
    }

    // Return whatever we have concluded.
    return c;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genome.GenomicRegion#toString()
   */
  @Override
  public String toString() {
    if (mText == null) {
      StringBuilder buffer = new StringBuilder();

      buffer.append(getType());
      buffer.append(" ").append(super.toString());
      buffer.append(" ").append(getStrand());
      buffer.append(" [");

      buffer.append(mPropertyMap.entrySet().stream().sorted(Map.Entry.comparingByKey())
          .map(entry -> entry.getKey() + "=" + entry.getValue()).collect(Collectors.joining(", ")));

      buffer.append("]");

      mText = buffer.toString();
    }

    return mText;
  }

  //
  // Static methods
  //

  /**
   * Returns the correct start based on the strand. Coordinates are stored on the
   * forward strand, so elements on the reverse strand must return their end as
   * their start.
   *
   * @param gene the gene
   * @return the genomic region
   */
  public static GenomicRegion tssRegion(GenomicRegion gene) {
    if (gene == null) {
      return null;
    }

    if (gene.mStrand == Strand.SENSE) {
      return new GenomicRegion(gene.mChr, gene.mStart, gene.mStart);
    } else {
      return new GenomicRegion(gene.mChr, gene.mEnd, gene.mEnd);
    }
  }

  /**
   * Tss dist.
   *
   * @param gene   the gene
   * @param region the region
   * @return the int
   */
  public static int tssDist(GenomicElement gene, GenomicRegion region) {
    GenomicRegion tssRegion = tssRegion(gene);

    if (gene.mStrand == Strand.SENSE) {
      return GenomicRegion.midDist(region, tssRegion);
    } else {
      return GenomicRegion.midDist(tssRegion, region);
    }
  }

  /**
   * Tss dist5p.
   *
   * @param gene   the gene
   * @param region the region
   * @return the int
   */
  public static int tssDist5p(GenomicElement gene, GenomicRegion region) {
    GenomicRegion tssRegion = tssRegion(gene);

    return GenomicRegion.midDist(region, tssRegion);
  }

  /**
   * To map.
   * 
   * @param <T>
   *
   * @param genes the genes
   * @return the iter map
   */
  public static <T extends GenomicElement> IterMap<Chromosome, Set<GenomicElement>> toMap(Collection<T> genes) {
    IterMap<Chromosome, Set<GenomicElement>> ret = DefaultTreeMap.create(new TreeSetCreator<GenomicElement>());

    for (GenomicElement g : genes) {
      ret.get(g.mChr).add(g);
    }

    return ret;
  }

  /**
   * Set the color of an element.
   * 
   * @param element
   * @param color
   * @param propogate
   */
  private static void setColor(GenomicElement element, Color color, boolean propogate) {
    Deque<GenomicElement> stack = new ArrayDeque<GenomicElement>();

    stack.push(element);

    GenomicElement e;

    while (!stack.isEmpty()) {
      e = stack.pop();
      e.mColor = color;

      if (propogate) {
        for (Entry<GenomicType, List<GenomicElement>> item : element.getChildren()) {
          for (GenomicElement c : item.getValue()) {
            stack.push(c);
          }
        }
      }
    }
  }

  /**
   * Get the distance from the mid point of a region to a gene accounting for the
   * strand.
   * 
   * @param gene
   * @param region
   * @return
   */
  public static int getTssMidDist(GenomicElement gene, GenomicRegion region) {
    int mid = GenomicRegion.mid(region);

    return getTssMidDist(gene, mid);
  }

  /**
   * Returns the distance of the mid to the gene tss. If the mid is downstream,
   * the value is positive.
   * 
   * @param gene
   * @param mid
   * @return
   */
  public static int getTssMidDist(GenomicElement gene, int mid) {
    if (gene.mStrand == Strand.SENSE) {
      return mid - gene.getStart();
    } else {
      return gene.getEnd() - mid;
    }
  }

  public Iterable<Entry<String, Object>> getProperties() {
    return mPropertyMap.entrySet().stream().sorted(Map.Entry.comparingByKey()).collect(Collectors.toList());
  }

}
