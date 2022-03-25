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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.jebtk.core.collections.IterMap;
import org.jebtk.core.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public abstract class GeneParser {
  public static final Logger LOG = LoggerFactory.getLogger(GeneParser.class);

  protected Set<GenomicType> mLevels = new HashSet<GenomicType>();

  /** Whether to add exons to gene structure */
  protected boolean mKeepExons = true;
  protected boolean mKeepTranscripts = true;
  protected Set<String> mMatchTags = new HashSet<String>();
  protected Set<String> mExcludeTags = new HashSet<String>();
  protected Set<String> mExcludeIds = new HashSet<String>();

  public GeneParser() {
    // setLevels(GeneType.GENE);
  }

  public GeneParser(GeneParser parser) {
    mLevels.addAll(parser.mLevels);
    mKeepExons = parser.mKeepExons;
    mMatchTags.addAll(parser.mMatchTags);
    mExcludeTags.addAll(parser.mExcludeTags);
  }

  public GeneParser setKeepExons(boolean keep) {
    GeneParser parser = create(this);

    parser._setKeepExons(keep);

    return parser;
  }

  protected void _setKeepExons(boolean keep) {
    mKeepExons = keep;
  }

  public GeneParser excludeIds(String id, String... ids) {
    GeneParser parser = create(this);

    parser.mExcludeIds.add(id);

    for (String i : ids) {
      parser.mExcludeIds.add(i);
    }

    return parser;
  }

  /**
   * Exclude entries matching given tags.
   * 
   * @param tag
   * @param tags
   * @return
   */
  public GeneParser excludeByTag(String tag, String... tags) {
    GeneParser parser = create(this);

    parser.mExcludeTags.add(tag);

    for (String t : tags) {
      parser.mExcludeTags.add(t);
    }

    return parser;
  }

  public GeneParser excludeByTag(Collection<String> excludeTags) {
    GeneParser parser = create(this);

    parser.mExcludeTags.addAll(excludeTags);

    return parser;
  }

  public GeneParser matchOnTag(String tag, String... tags) {
    GeneParser parser = create(this);

    parser.mMatchTags.add(tag);

    for (String t : tags) {
      parser.mMatchTags.add(t);
    }

    return parser;
  }

  public GeneParser setLevels(GenomicType level, GenomicType... levels) {
    GeneParser parser = create(this);

    parser._setLevels(level, levels);

    return parser;
  }

  protected void _setLevels(GenomicType level, GenomicType... levels) {
    mLevels.clear();
    mLevels.add(level);

    for (GenomicType l : levels) {
      mLevels.add(l);
    }
  }

  public GeneParser setLevels(Collection<GenomicType> levels) {
    GeneParser parser = create(this);

    parser._setLevels(levels);

    return parser;
  }

  protected void _setLevels(Collection<GenomicType> levels) {
    mLevels.clear();

    _addLevels(levels);
  }

  public GeneParser addLevels(Collection<GenomicType> levels) {
    GeneParser parser = create(this);

    parser._addLevels(levels);

    return parser;
  }

  protected void _addLevels(Collection<GenomicType> levels) {
    mLevels.addAll(levels);
  }

  public abstract GeneParser create(GeneParser parser);

  public GenesDB parse(Path file, final Genome genome) throws IOException {
    GenesDB genes = new FixedGapGenes(genome);

    parse(file, genome, genes);

    return genes;
  }

  protected GenesDB parse(Path file, Genome genome, BufferedReader reader) throws IOException {
    GenesDB genes = new FixedGapGenes(genome);

    parse(file, reader, genome, genes);

    return genes;
  }

  public void parse(Path file, final Genome genome, GenesDB genes) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    try {
      parse(file, reader, genome, genes);
    } finally {
      reader.close();
    }
  }

  protected void parse(Path file, BufferedReader reader, final Genome genome, GenesDB genes) throws IOException {
    // Do nothing
  }

  public void parse(Path file, final Genome genome, Chromosome chr, GenesDB genes) throws IOException {
    parse(file, genome, genes);
  }

  public GenesDB parse(Path file, final Genome genome, int window) throws IOException {
    GenesDB genes = new FixedGapGenes(genome);

    parse(file, genes, genome, window);

    return genes;
  }

  protected GenesDB parse(Path file, BufferedReader reader, final Genome genome, int window) throws IOException {
    GenesDB genes = new FixedGapGenes(genome);

    parse(file, reader, genes, genome, window);

    return genes;
  }

  public void parse(Path file, GenesDB genes, final Genome genome, int window) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    try {
      parse(file, reader, genes, genome, window);
    } finally {
      reader.close();
    }
  }

  protected void parse(Path file, BufferedReader reader, GenesDB genes, final Genome genome, int window)
      throws IOException {
    parse(file, reader, genome, genes);
  }

  public Map<String, Set<String>> idMap(Path file, String id1, String id2) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    Map<String, Set<String>> ret;

    try {
      ret = idMap(file, reader, id1, id2);
    } finally {
      reader.close();
    }

    return ret;
  }

  /**
   * Parses the gene table.
   *
   * @param reader the reader
   * @return the genes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public abstract Map<String, Set<String>> idMap(Path file, BufferedReader reader, String id1, String id2)
      throws IOException;

  public boolean containsLevel(GenomicType level) {
    if (mLevels.size() == 0) {
      return true;
    }

    return mLevels.contains(level);
  }

  public static GenomicEntity addAttributes(GenomicType type, final GenomicRegion region,
      final IterMap<String, String> attributeMap) {

    GenomicEntity gene = new GenomicEntity(type, region);

    // Add the ids
    for (Entry<String, String> f : attributeMap) {
      gene.setProperty(f.getKey(), f.getValue());
    }

    // If there are any tags

    // genes.add(gene, gene);

    return gene;
  }
}