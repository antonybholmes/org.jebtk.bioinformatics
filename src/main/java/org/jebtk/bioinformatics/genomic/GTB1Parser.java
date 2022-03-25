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
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.TreeSetCreator;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public class GTB1Parser extends GTBParser {

  public GTB1Parser() {
    // _setLevels(GeneType.GENE);
  }

  public GTB1Parser(GeneParser parser) {
    super(parser);
  }

  /**
   * Parses the gene table.
   *
   * @param reader the reader
   * @return the genes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  @Override
  protected void parse(Path file, BufferedReader reader, final Genome genome, GenesDB genes) throws IOException {
    LOG.info("Parsing GTB file {}, levels: {}...", file, mLevels);

    String line;
    List<String> tokens;

    GenomicEntity gene = null;

    // Skip header
    reader.readLine();

    Splitter splitter = Splitter.on(';');

    // Add the exons
    boolean hasExonLevel = containsLevel(GenomicType.EXON);

    while ((line = reader.readLine()) != null) {
      if (Io.isEmptyLine(line)) {
        continue;
      }

      boolean add = true;

      tokens = Splitter.onTab().text(line);

      Chromosome chr = ChromosomeService.getInstance().chr(genome, tokens.get(0));

      // Skip random and unofficial chromosomes
      if (chr.toString().contains("_")) {
        continue;
      }

      Strand strand = Strand.parse(tokens.get(1));
      int start = Integer.parseInt(tokens.get(2));
      int end = Integer.parseInt(tokens.get(3));

      // int exonCount = Integer.parseInt(tokens.get(4));

      // Because of the UCSC using zero based start and one
      // based end, we need to increment the start by 1

      List<Integer> starts = TextUtils.splitInts(tokens.get(5), TextUtils.SEMI_COLON_DELIMITER);

      List<Integer> ends = TextUtils.splitInts(tokens.get(6), TextUtils.SEMI_COLON_DELIMITER);

      List<String> tags = null;

      if (tokens.size() > 8) {
        tags = TextUtils.removeNA(splitter.text(tokens.get(8)));

        if (mExcludeTags.size() > 0) {
          for (String tag : tags) {
            if (mExcludeTags.contains(tag)) {
              add = false;
              break;
            }
          }
        }

        if (mMatchTags.size() > 0) {
          add = false;

          for (String tag : tags) {
            if (mMatchTags.contains(tag)) {
              add = true;
              break;
            }
          }
        }
      }

      if (!add) {
        continue;
      }

      IterMap<String, String> attributeMap = getAttributes(splitter, tokens.get(7));

      // Create the gene

      gene = addAttributes(GenomicType.TRANSCRIPT, GenomicRegion.create(chr, start, end, strand), attributeMap);

      if (containsLevel(GenomicType.TRANSCRIPT)) {
        genes.add(gene);
      }

      if (hasExonLevel || mKeepExons) {
        for (int i = 0; i < starts.size(); ++i) {
          // Again correct for the ucsc
          GenomicRegion region = GenomicRegion.create(chr, starts.get(i) + 1, ends.get(i), strand);

          GenomicEntity exon = addAttributes(GenomicType.EXON, region, attributeMap);

          if (mKeepExons) {
            if (gene != null) {
              gene.addChild(exon);
            }
          }

          if (hasExonLevel) {

            exon.setParent(gene);

            genes.add(exon);
          }
        }
      }
    }
  }

  @Override
  public Map<String, Set<String>> idMap(Path file, BufferedReader reader, String id1, String id2) throws IOException {
    LOG.info("Creating id map from GTB file {}, levels: {}...", file, mLevels);

    Map<String, Set<String>> ret = DefaultTreeMap.create(new TreeSetCreator<String>());

    String line;
    List<String> tokens;

    // Skip header
    reader.readLine();

    Splitter splitter = Splitter.on(';');

    while ((line = reader.readLine()) != null) {
      if (Io.isEmptyLine(line)) {
        continue;
      }

      boolean add = true;

      tokens = Splitter.onTab().text(line);

      Chromosome chr = ChromosomeService.getInstance().guessChr(file, tokens.get(0));

      // Skip random and unofficial chromosomes
      if (chr.toString().contains("_")) {
        continue;
      }

      if (tokens.size() > 8) {
        List<String> tags = TextUtils.removeNA(splitter.text(tokens.get(8)));

        if (mExcludeTags.size() > 0) {
          for (String tag : tags) {
            if (mExcludeTags.contains(tag)) {
              add = false;
              break;
            }
          }
        }

        if (mMatchTags.size() > 0) {
          add = false;

          for (String tag : tags) {
            if (mMatchTags.contains(tag)) {
              add = true;
              break;
            }
          }
        }
      }

      if (!add) {
        continue;
      }

      IterMap<String, String> attributeMap = getAttributes(splitter, tokens.get(7));

      String name1 = attributeMap.get(id1);
      String name2 = attributeMap.get(id2);

      // System.err.println("id " + id1 + " " + name1 + " " + id2 + " " +
      // name2);

      if (name1 != null && name2 != null) {
        ret.get(name1).add(name2);
      }
    }

    return ret;
  }

  @Override
  public GeneParser create(GeneParser parser) {
    return new GTB1Parser(parser);
  }

}