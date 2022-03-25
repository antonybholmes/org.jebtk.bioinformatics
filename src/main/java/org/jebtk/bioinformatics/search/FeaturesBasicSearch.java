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

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import org.jebtk.bioinformatics.Search;
import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.ChromosomeService;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * The class FeaturesBasicSearch.
 */
public class FeaturesBasicSearch extends AbstractFeaturesSearch {
  // protected Map<Short, List<Feature>> allLocations =
  // new HashMap<Short, List<Feature>>();

  /**
   * The all locations.
   */
  protected List<List<Feature>> allLocations = new CopyOnWriteArrayList<List<Feature>>();

  /**
   * The member file.
   */
  protected Path mFile;

  /**
   * The size.
   */
  private int size = 0;

  /**
   * Instantiates a new features basic search.
   *
   * @param name        the name
   * @param description the description
   * @param file        the file
   */
  public FeaturesBasicSearch(String name, String description, Path file) {
    super(name, description);

    mFile = file;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * cacheFeatures()
   */
  public final void cacheFeatures() {
    allLocations.clear();

    for (int i = 0; i < 26; ++i) {
      allLocations.add(new ArrayList<Feature>());
    }

    try {
      BufferedReader reader = FileUtils.newBufferedReader(mFile);

      String line;

      try {
        // skip header
        line = reader.readLine();

        while ((line = reader.readLine()) != null) {
          if (Io.isEmptyLine(line)) {
            continue;
          }

          List<String> row = TextUtils.fastSplit(line, TextUtils.TAB_DELIMITER);

          Feature feature = new Feature(row.get(0), ChromosomeService.getInstance().guessChr(mFile, row.get(1)),
              Integer.parseInt(row.get(2)), Integer.parseInt(row.get(3)));
          // feature.type = type;

          // System.err.println(line + " chr:" + feature.getChromosome());

          allLocations.get(feature.getChr().getId()).add(feature);
          mFeatureByName.put(feature.getName(), feature);

          ++size;
        }
      } finally {
        reader.close();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome)
   */
  public final List<Feature> getFeatures(Chromosome chromosome) {
    if (allLocations.size() == 0) {
      cacheFeatures();
    }

    return allLocations.get(chromosome.getId());
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * freeMemory( )
   */
  public void freeMemory() {
    super.freeMemory();

    allLocations.clear();
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome, int, int)
   */
  @Override
  public List<Feature> getFeatures(Chromosome chromosome, int start, int endLocation) {
    if (allLocations.size() == 0) {
      cacheFeatures();
    }

    // System.err.println(getName() + " " + chromosome + ":loc:" + startLocation
    // + "
    // " + endLocation);

    // go through the table

    int startIndex;
    int endIndex;

    List<Feature> features = new ArrayList<Feature>();

    List<Feature> locations = allLocations.get(chromosome.getId());

    if (locations.size() == 0) {
      return features;
    }
    // endLocations = allEndLocations.get(chromosome);

    int last = locations.size() - 1;

    if (start > locations.get(last).getEnd() || endLocation < locations.get(0).getStart()) {
      System.out.println(start + ":" + locations.get(0).getStart() + ":" + endLocation + ":"
          + locations.get(locations.size() - 1).getEnd());

      // the range is clearly not within the feature set so don't even bother to
      // look
      return features;
    }

    if (start <= locations.get(0).getStart()) {
      startIndex = 0;
    } else {
      startIndex = Search.findInner(start, locations);
    }

    if (startIndex == -1) {
      // there was no valid start within the boundaries so no point
      // continuing

      return features;
    }

    if (endLocation >= locations.get(last).getStart()) {
      endIndex = last;
    } else {
      endIndex = Search.findOuter(endLocation, locations);
    }

    if (endIndex == -1) {
      // there was no valid start within the boundaries so no point
      // continuing

      return features;
    }

    for (int i = startIndex; i <= endIndex; ++i) {
      features.add(locations.get(i));
    }

    return features;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome, int)
   */
  @Override
  public List<Feature> getFeatures(Chromosome chromosome, int location) {
    if (allLocations.size() == 0) {
      cacheFeatures();
    }

    List<Feature> features = new ArrayList<Feature>();

    List<Feature> locations = allLocations.get(chromosome.getId());

    if (locations.size() == 0) {
      return features;
    }

    for (int i = 0; i < locations.size(); ++i) {
      if (features.get(i).getStart() <= location && features.get(i).getEnd() >= location) {
        features.add(features.get(i));
      }
    }

    return features;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#size()
   */
  @Override
  public int size() {
    return size;
  }
}
