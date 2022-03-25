/**
 * Copyright 2017 Antony Holmes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.jebtk.bioinformatics.genomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.bioinformatics.gapsearch.FixedGapSearch;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.collections.TreeSetCreator;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

/**
 * The Class GFF.
 */
public class GFF3Parser extends GeneParser {

  public static final String LEVEL_GENE = "gene";
  public static final String LEVEL_TRANSCRIPT = "transcript";
  public static final String LEVEL_EXON = "exon";

  public GFF3Parser() {
    // _setKeepExons(false);
  }

  public GFF3Parser(GeneParser parser) {
    super(parser);
  }

  @Override
  protected void parse(Path file, BufferedReader reader, final Genome genome, GenesDB genes) throws IOException {
    String line;
    List<String> tokens;
    GenomicType type;
    int start;
    int end;
    Strand strand;
    IterMap<String, String> attributeMap;

    Splitter splitter = Splitter.onTab();

    GenomicEntity gene = null;
    GenomicEntity transcript = null;

    try {
      while ((line = reader.readLine()) != null) {
        // System.err.println(line);

        tokens = splitter.text(line);

        Chromosome chr = ChromosomeService.getInstance().chr(genome, tokens.get(0));

        type = GenomicType.parse(tokens.get(2));
        start = Integer.parseInt(tokens.get(3));
        end = Integer.parseInt(tokens.get(4));
        strand = Strand.parse(tokens.get(6));

        // Must be at least the current level or lower
        if (!containsLevel(type)) {
          continue;
        }

        attributeMap = parseAttributes(tokens.get(8));

        GenomicRegion region = GenomicRegion.create(chr, start, end, strand);

        switch (type) {
        case GENE:
          gene = addAttributes(GenomicType.GENE, region, attributeMap);

          genes.add(gene);

          break;
        case TRANSCRIPT:
          // Gene
          // gene = new Gene(region).setSymbol(name);

          transcript = addAttributes(GenomicType.TRANSCRIPT, region, attributeMap);

          if (gene != null) {
            if (mKeepTranscripts) {
              gene.addChild(transcript);
            }

            transcript.setParent(gene);
          }

          genes.add(transcript);

          break;
        case EXON:
          GenomicEntity exon = addAttributes(GenomicType.EXON, region, attributeMap);

          if (transcript != null) {
            // Add to the current gene

            if (mKeepExons) {
              transcript.addChild(exon);
            }

            transcript.setParent(gene);
          }

          genes.add(exon);

          break;
        default:
          break;
        }
      }
    } finally {
      reader.close();
    }
  }

  @Override
  public Map<String, Set<String>> idMap(Path file, BufferedReader reader, String id1, String id2) throws IOException {
    Map<String, Set<String>> ret = DefaultTreeMap.create(new TreeSetCreator<String>());

    String line;
    List<String> tokens;

    GenomicType type;

    Map<String, String> attributeMap;

    Splitter splitter = Splitter.onTab();

    try {
      while ((line = reader.readLine()) != null) {
        // System.err.println(line);

        tokens = splitter.text(line);

        type = GenomicType.parse(tokens.get(2));

        // Must be at least the current level or lower
        if (!containsLevel(type)) {
          continue;
        }

        attributeMap = parseAttributes(tokens.get(8));

        String name1 = attributeMap.get(id1);
        String name2 = attributeMap.get(id2);

        if (name1 != null && name2 != null) {
          ret.get(name1).add(name2);
        }
      }
    } finally {
      reader.close();
    }

    return ret;
  }

  /**
   * Parses the attributes.
   *
   * @param attributes the attributes
   * @return the map
   */
  public static IterMap<String, String> parseAttributes(String attributes) {
    // System.err.println("attributes " + attributes);

    IterMap<String, String> ret = new IterTreeMap<String, String>();

    List<String> tokens = Splitter.on(TextUtils.SEMI_COLON_DELIMITER_CHAR).text(formatAttributes(attributes));

    for (String token : tokens) {
      // System.err.println("gff3 " + token);
      List<String> values = Splitter.on('=').text(token);

      ret.put(values.get(0), values.get(1));
    }

    return ret;
  }

  public static String formatAttributes(String attributes) {
    String ret = attributes.trim().replace("  ", TextUtils.SPACE_DELIMITER)
        .replace(" ;", TextUtils.SEMI_COLON_DELIMITER).replace("; ", TextUtils.SEMI_COLON_DELIMITER)
        .replace("\";", TextUtils.SEMI_COLON_DELIMITER).replace(" \"", "=").replace("\"", TextUtils.EMPTY_STRING);

    // Remove last semi-colon if there is one
    if (ret.charAt(ret.length() - 1) == TextUtils.SEMI_COLON_DELIMITER_CHAR) {
      ret = ret.substring(0, ret.length() - 1);
    }

    return ret;
  }

  /**
   * Parses the GFF 3 attributes.
   *
   * @param tokens the tokens
   * @return the map
   */
  public static Map<String, String> parseGFF3Attributes(List<String> tokens) {
    return parseGFF3Attributes(tokens.get(8));
  }

  /**
   * Parse gff3 formatted attributes, e.g. gene_id="BCL6";exon="2";
   *
   * @param attributes the attributes
   * @return the map
   */
  public static Map<String, String> parseGFF3Attributes(String attributes) {
    // System.err.println("attributes " + attributes);

    Map<String, String> ret = new HashMap<String, String>();

    List<String> tokens = Splitter.on(';').text(attributes);

    Splitter equalsSplitter = Splitter.on('=');

    for (String token : tokens) {
      List<String> values = equalsSplitter.text(token);

      if (values.size() > 1) {
        String name = values.get(0).trim();
        String value = values.get(1).trim().replace("\"", "");
        ret.put(name, value);

        // System.err.println("attribute " + name + " " + value);
      }
    }

    return ret;
  }

  /**
   * Returns an attribute name where underscores are replaced with spaces and the
   * name is converted to sentence case to make it more presentable in a table
   * header etc.
   *
   * @param attribute the attribute
   * @return the string
   */
  public static String formatAttributeName(String attribute) {
    return TextUtils.titleCase(attribute.replace("_", " "));
  }

  public static FixedGapSearch<Gene> GFFToGapSearch(List<Gene> features) {
    FixedGapSearch<Gene> ret = new FixedGapSearch<Gene>(1000);

    for (Gene gene : features) {
      ret.add(gene, gene);
    }

    return ret;
  }

  @Override
  public GeneParser create(GeneParser parser) {
    return new GFF3Parser(parser);
  }

  public static List<String> gff3IdTypes(Path gff3, GenomicType type) throws IOException {

    BufferedReader reader = FileUtils.newBufferedReader(gff3);

    List<String> ret = new ArrayList<String>(20);

    String line;

    List<String> tokens;

    Splitter splitter = Splitter.onTab();

    GenomicType t;

    try {
      while ((line = reader.readLine()) != null) {
        tokens = splitter.text(line);

        t = GenomicType.parse(tokens.get(2));

        if (t == type) {
          Map<String, String> attributes = parseGFF3Attributes(TextUtils.tabSplit(line));

          for (String attribute : CollectionUtils.sortKeys(attributes)) {

            // Skip some predefined fields that will not be of interest to the
            // user
            if (attribute.contains("Parent")) {
              continue;
            }

            ret.add(attribute); // formatAttributeName(attribute));
          }

          break;
        }
      }
    } finally {
      reader.close();
    }

    Collections.sort(ret);

    return ret;
  }
}
