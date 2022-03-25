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
package org.jebtk.bioinformatics.search;

import java.nio.file.Path;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.Search;
import org.jebtk.bioinformatics.genomic.Chromosome;

/**
 * Provides binary search within an ordered list of features with the proviso
 * that none of the features have starting or stopping overlaps.
 *
 * @author Antony Holmes
 *
 */
public class FeaturesBinarySearch extends FeaturesBasicSearch {

  /**
   * Instantiates a new features binary search.
   *
   * @param name        the name
   * @param description the description
   * @param file        the file
   */
  public FeaturesBinarySearch(String name, String description, Path file) {
    super(name, description, file);
  }

  /**
   * Gets the features.
   *
   * @param startLocation the start location
   * @param endLocation   the end location
   * @param chromosome    the chromosome
   * @return the features
   * @throws ParseException the parse exception
   */
  public List<Feature> getFeatures(int startLocation, int endLocation, Chromosome chromosome) {
    // System.out.println("loc:" + startLocation + " " + endLocation);

    // go through the table

    int startIndex;
    int endIndex;

    List<Feature> features = new ArrayList<Feature>(100);

    List<Feature> locations = allLocations.get(chromosome.getId());

    if (locations == null) {
      // System.out.println("ropey chromosome " + chromosome);
      return features;
    }
    // endLocations = allEndLocations.get(chromosome);

    // System.out.println(startLocation + ":" + locations.get(0).getStart() +
    // ":" +
    // endLocation + ":" + locations.get(locations.size() - 1).getEnd());

    int last = locations.size() - 1;

    if (startLocation > locations.get(last).getEnd() || endLocation < locations.get(0).getStart()) {
      // the range is clearly not within the feature set so don't even bother to
      // look
      return features;
    }

    // the range is within some part of the features

    // startIndex = Search.binaryRangeSearch(startLocation, locations,
    // SearchType.STRICT);

    if (startLocation <= locations.get(0).getStart()) {
      startIndex = 0;
    } else {
      startIndex = Search.binaryRangeSearchInner(startLocation, locations);
    }

    // if (startIndex == -1) {
    // means there was nothing in the range
    // return features;
    // }

    if (endLocation >= locations.get(last).getStart()) {
      endIndex = last;
    } else {
      endIndex = Search.binaryRangeSearchOuter(endLocation, locations);
    }

    // if (endIndex == -1) {
    // means there was nothing in the range
    // return features;
    // }

    // System.out.println(startIndex + " " + endIndex);

    for (int i = startIndex; i <= endIndex; ++i) {
      features.add(locations.get(i));
    }

    return features;
  }

  /**
   * The main method.
   *
   * @param args the arguments
   */
  public static final void main(String[] args) {
    // convertToBin(new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/probes.txt"),
    // new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/probes.bin"),
    // Ui.SMALL_ICON_SIZE);

    // convertToBin(new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/cytobands_condensed.txt"),
    // new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/cytobands_condensed.bin"),
    // 10);

    // convertToBin(new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/cnv.txt"),
    // new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/cnv.bin"),
    // 2);

    /*
     * convertToBin(new File(
     * "/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/genes.txt"
     * ), new File(
     * "/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/genes.bin"
     * ), Ui.SMALL_ICON_SIZE);
     * 
     * convertToBin(new File(
     * "/home/antony/university/columbia/rdf/copy_number/resources/mirbase/annotations/17/mirs.txt"
     * ), new File(
     * "/home/antony/university/columbia/rdf/copy_number/resources/mirbase/annotations/17/mirs.bin"
     * ), 20);
     */

    // convertToBin(new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/mirs.txt"),
    // new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/mirs.bin"),
    // 20);

    // convertToBin(new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/genes_ccds.txt"),
    // new
    // File("/home/antony/university/columbia/rdf/copy_number/resources/ucsc/references/assembly/hg18/genes_ccds.bin"),
    // Ui.SMALL_ICON_SIZE);
    /*
     * Features features = new Features("probes");
     * 
     * features.cacheFeaturesBin(new File(
     * "/home/antony/university/columbia/rdf/copy_number/resources/affymetrix/references/snp_6/annotations/30/probes.bin"
     * ));
     * 
     * List<Feature> f = features.getFeatures(1, 1000000, 1);
     * 
     * for (Feature l : f) { System.out.println(l.feature + " " + l.group + " " +
     * l.start); }
     */
  }

  /*
   * public static final void convertToBin(File in, File out, int stringSize) {
   * String line;
   * 
   * int size = 0;
   * 
   * BufferedReader reader = null;
   * 
   * Map<Short, Integer> groupCounts = new HashMap<Short, Integer>();
   * 
   * short chromosome;
   * 
   * try { reader = new BufferedReader(new FileReader(in));
   * 
   * // skip header line = reader.readLine();
   * 
   * try { while ((line = reader.readLine()) != null) { if (Io.isEmptyLine(line))
   * { continue; }
   * 
   * List<String> row = Text.fastSplitRemoveQuotes(line);
   * 
   * chromosome = Chromosome.parse(row.get(1));
   * 
   * if (!groupCounts.containsKey(chromosome)) { groupCounts.put(chromosome, 0); }
   * 
   * groupCounts.put(chromosome, groupCounts.get(chromosome) + 1);
   * 
   * ++size; } } finally { reader.close(); } } catch (Exception e) {
   * e.printStackTrace(); }
   * 
   * // now actually parse the file try { reader = new BufferedReader(new
   * FileReader(in));
   * 
   * DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new
   * FileOutputStream(out)));
   * 
   * try { //DataOutputStream dos2;
   * 
   * // first write out the string size so we know how to skip dos.writeInt(size);
   * 
   * // first write out the string size so we know how to skip
   * dos.write(stringSize);
   * 
   * String name;
   * 
   * int startLocation; int endLocation;
   * 
   * // skip header line = reader.readLine();
   * 
   * while ((line = reader.readLine()) != null) { if (Io.isEmptyLine(line)) {
   * continue; }
   * 
   * List<String> row = Text.fastSplitRemoveQuotes(line);
   * 
   * name = row.get(0);
   * 
   * chromosome = Chromosome.parse(row.get(1));
   * 
   * startLocation = Integer.parseInt(row.get(2));
   * 
   * endLocation = Integer.parseInt(row.get(3));
   * 
   * // write out feature in padded byte array
   * 
   * dos.write(Arrays.copyOf(name.getBytes(), stringSize)); dos.write(chromosome);
   * dos.writeInt(startLocation); dos.writeInt(endLocation); } } finally {
   * reader.close(); dos.close(); } } catch (Exception e) { e.printStackTrace(); }
   * }
   */
}
