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
package org.jebtk.bioinformatics.pathway;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.settings.SettingsService;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.statistics.Hypergeometric;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The class Pathway.
 */
public class Pathway {

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(Pathway.class);

  // Match genepattern
  private static final int TOTAL_GENES = SettingsService.getInstance().getInt("bioinformatics.pathway.total-genes");

  /**
   * Analysis.
   *
   * @param names             the names
   * @param collections       the collections
   * @param refseqConversion  the refseq conversion
   * @param ensemblConversion the ensembl conversion
   * @param maxFdr            the max fdr
   * @param out               the out
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void analysis(List<String> names, Collection<GeneSet> collections, IdToSymbol refseqConversion,
      IdToSymbol ensemblConversion, double maxFdr, Path out) throws IOException {
    analysis(names, collections, refseqConversion, ensemblConversion, maxFdr, out, null);
  }

  /**
   * Analysis.
   *
   * @param names             the names
   * @param collections       the collections
   * @param refseqConversion  the conversion
   * @param ensemblConversion the ensembl conversion
   * @param maxFdr            the max fdr
   * @param out               the out
   * @param tableFile         the table file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void analysis(List<String> names, Collection<GeneSet> collections, IdToSymbol refseqConversion,
      IdToSymbol ensemblConversion, double maxFdr, Path out, Path tableFile) throws IOException {

    LOG.info("Max FDR = {}", maxFdr);

    int totalGeneSets = collections.size();

    LOG.info("Gene set collections = {}", collections.size());

    LOG.info("Total gene sets = {}", totalGeneSets);

    // determine which reference we are using

    int matches = 0;

    for (int i = 0; i < Math.min(10, names.size()); ++i) {
      if (refseqConversion.convert(names.get(i)) != null) {
        matches += 1;
      }
    }

    IdToSymbol conversion = null;

    if (matches > 2) {
      conversion = refseqConversion;
    } else {
      matches = 0;

      for (int i = 0; i < Math.min(5, names.size()); ++i) {
        if (ensemblConversion.convert(names.get(i)) != null) {
          matches += 1;
        }
      }

      if (matches > 2) {
        conversion = ensemblConversion;
      }
    }

    if (conversion == null) {
      return;
    }

    // See if the symbols are valid

    List<String> validSymbols = new ArrayList<String>();

    for (String id : names) {
      String symbol = conversion.convert(id);

      if (symbol != null) {
        validSymbols.add(symbol);
      }
    }

    int totalGenes = TOTAL_GENES; // conversion.getSymbolCount();

    LOG.info("Gene universe = {}", totalGenes);

    // Lets do some overlaps

    Hypergeometric hyperg = new Hypergeometric();

    Map<Double, Set<GeneSet>> pvalues = new HashMap<Double, Set<GeneSet>>();

    Map<GeneSet, Double> geneSetPMap = new HashMap<GeneSet, Double>();

    Map<GeneSet, Integer> overlaps = new HashMap<GeneSet, Integer>();

    Map<GeneSet, Double> overlapRatios = new HashMap<GeneSet, Double>();

    for (GeneSet geneSet : collections) {
      int overlap = 0;

      for (String symbol : validSymbols) {
        if (geneSet.contains(symbol)) {
          ++overlap;
        }
      }

      overlaps.put(geneSet, overlap);

      double overlapGeneSetRatio = overlap / (double) geneSet.size();
      overlapRatios.put(geneSet, overlapGeneSetRatio);

      /*
       * double p = hyperg.pdf(overlap, validSymbols.size(), geneSet.size(),
       * totalGenes);
       */

      // Default to assume the enrichment is not significant.
      double p = 1.0;

      // Sum of probabilities of seeing overlap
      // or more genes by chance, If we use the
      // EASE method, we also reduce the overlap
      // by one.
      // https://david.ncifcrf.gov/helps/functional_annotation.html

      if (overlap > 1) {
        p = 1 - hyperg.cdf(overlap - 2, geneSet.size(), validSymbols.size(), totalGenes);

        p = Mathematics.bound(p, 0, 1);
      }

      if (!pvalues.containsKey(p)) {
        pvalues.put(p, new TreeSet<GeneSet>());
      }

      pvalues.get(p).add(geneSet);

      geneSetPMap.put(geneSet, p);
    }

    List<Double> sortedPValues = CollectionUtils.sortKeys(pvalues);

    Map<Double, Double> fdrs = new HashMap<Double, Double>();

    for (int i = 0; i < sortedPValues.size(); ++i) {
      double p = sortedPValues.get(i);

      int rank = i + 1;

      double fdr = Mathematics.bound(p * totalGeneSets / rank, 0, 1);

      fdrs.put(fdr, p);
    }

    // Find the genesets with acceptable q-values

    List<GeneSet> acceptedGeneSets = new ArrayList<GeneSet>();
    List<Double> acceptedFdrs = new ArrayList<Double>();

    for (double fdr : CollectionUtils.sort(fdrs.keySet())) {
      double p = fdrs.get(fdr);

      System.err.println("fdr " + fdr + " " + maxFdr);

      if (fdr > maxFdr) {
        continue;
      }

      acceptedFdrs.add(fdr);

      for (GeneSet geneSet : CollectionUtils.sort(pvalues.get(p))) {
        acceptedGeneSets.add(geneSet);
      }
    }

    Map<String, List<Boolean>> overlapMatrix = new HashMap<String, List<Boolean>>();

    for (String gene : validSymbols) {
      List<Boolean> in = new ArrayList<Boolean>();

      for (GeneSet geneSet : acceptedGeneSets) {
        in.add(geneSet.contains(gene));
      }

      overlapMatrix.put(gene, in);
    }

    //
    // Write out the results
    //

    BufferedWriter writer = FileUtils.newBufferedWriter(out);

    // First write the collection names

    DecimalFormat decimalFormat = new DecimalFormat("0.00");
    DecimalFormat scientificFormat = new DecimalFormat("0.##E0");

    try {
      Set<String> collectionNames = new TreeSet<String>();

      for (GeneSet geneSet : collections) {
        collectionNames.add(geneSet.getCollectionName());
      }

      writer.write("collections");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(TextUtils.scJoin(CollectionUtils.toList(collectionNames)));
      writer.newLine();
      writer.write("gene_sets");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(Integer.toString(totalGeneSets));
      writer.newLine();
      writer.write("valid_gene_symbols");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(Integer.toString(validSymbols.size()));
      writer.newLine();
      writer.write("genes");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(Integer.toString(totalGenes));
      writer.newLine();
      writer.newLine();

      writer.write("gene_set_name");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write("genes_in_gene_set");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write("genes_in_overlap");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write("overlap_gene_set_ratio");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write("p-value");
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write("fdr");
      writer.newLine();

      for (double fdr : acceptedFdrs) {
        double p = fdrs.get(fdr);

        for (GeneSet geneSet : CollectionUtils.sort(pvalues.get(p))) {
          int overlap = overlaps.get(geneSet);
          double overlapRatio = overlapRatios.get(geneSet);

          // System.err.println(geneSet.getName() + " " + overlap + " " +
          // validSymbols.size() + " " + geneSet.size() + " " + totalGenes + " "
          // + p);

          writer.write(geneSet.getName());
          writer.write(TextUtils.TAB_DELIMITER);
          writer.write(Integer.toString(geneSet.size()));
          writer.write(TextUtils.TAB_DELIMITER);
          writer.write(Integer.toString(overlap));
          writer.write(TextUtils.TAB_DELIMITER);
          // writer.write(Double.toString(overlapRatio));
          writer.write(decimalFormat.format(overlapRatio));
          writer.write(TextUtils.TAB_DELIMITER);
          // writer.write(Double.toString(p));
          writer.write(scientificFormat.format(p));
          writer.write(TextUtils.TAB_DELIMITER);
          // writer.write(Double.toString(fdr));
          writer.write(scientificFormat.format(fdr));

          writer.newLine();
        }
      }

      writer.newLine();

      writer.write("Gene Symbol");

      for (GeneSet geneSet : acceptedGeneSets) {
        writer.write("\t");
        writer.write(geneSet.getName());
      }

      writer.newLine();

      for (String gene : CollectionUtils.sort(overlapMatrix.keySet())) {
        // writer.write(gene);
        // writer.write(TextUtils.TAB_DELIMITER);
        writer.write(gene.toUpperCase());

        for (boolean in : overlapMatrix.get(gene)) {
          writer.write("\t");

          if (in) {
            writer.write("1");
          } else {
            writer.write("0");
          }
        }

        writer.newLine();
      }

    } finally {
      writer.close();
    }

    if (tableFile == null) {
      return;
    }

    //
    // Create file suitable for heatmap
    //

    writer = FileUtils.newBufferedWriter(tableFile);

    try {
      writer.write("Gene Symbol\tDescription");

      for (GeneSet geneSet : acceptedGeneSets) {
        writer.write("\t");
        writer.write(geneSet.getName());
      }

      writer.newLine();

      for (String gene : CollectionUtils.sort(overlapMatrix.keySet())) {
        writer.write(gene.toUpperCase());
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(gene.toUpperCase());

        for (boolean in : overlapMatrix.get(gene)) {
          writer.write("\t");

          if (in) {
            writer.write("3");
          } else {
            writer.write("-3");
          }
        }

        writer.newLine();
      }
    } finally {
      writer.close();
    }
  }
}
