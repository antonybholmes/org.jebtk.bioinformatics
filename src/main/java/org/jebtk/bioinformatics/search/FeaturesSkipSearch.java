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
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.Search;
import org.jebtk.bioinformatics.genomic.Chromosome;

/**
 * The class FeaturesSkipSearch.
 */
public class FeaturesSkipSearch extends FeaturesBinarySearch {

  /**
   * The skip.
   */
  private static int SKIP = 1000;

  /**
   * Instantiates a new features skip search.
   *
   * @param name        the name
   * @param description the description
   * @param file        the file
   */
  public FeaturesSkipSearch(String name, String description, Path file) {
    super(name, description, file);
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.FeaturesBinarySearch#
   * getFeatures( int, int, short)
   */
  @Override
  public final List<Feature> getFeatures(int startLocation, int endLocation, Chromosome chromosome) {
    // System.out.println(chromosome + ":loc:" + startLocation + " " +
    // endLocation);

    // go through the table

    int startIndex;
    int endIndex;

    List<Feature> features = new ArrayList<Feature>();

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
    } else if (startLocation <= locations.get(0).getStart() && endLocation >= locations.get(last).getStart()) {
      // the range spans the entire set of features

      startIndex = 0;
      endIndex = last;

      // System.out.println("f:" + last + " " + locations.size() + " " +
      // startIndex +
      // " " + endIndex + " " + locations.get(startIndex).getName() + " " +
      // locations.get(endIndex).getName());
    } else {
      // the range is within some part of the features

      // System.out.println(startLocation + " " + locations.get(0).getStart() +
      // " " +
      // locations.get(0).getEnd() + " " + endLocation + " " +
      // locations.get(last).getStart() + " " + locations.get(last).getEnd());

      if (startLocation <= locations.get(0).getStart()) {
        startIndex = 0;
      } else {
        startIndex = Search.jumpSearchInner(startLocation, locations, SKIP);
      }

      if (startIndex == -1) {
        // there was no valid start within the boundaries so no point
        // continuing
        return features;
      }

      if (endLocation >= locations.get(last).getStart()) {
        endIndex = last;
      } else {
        endIndex = Search.jumpSearchOuter(endLocation, locations, SKIP);
      }

      if (endIndex == -1) {
        // there was no valid start within the boundaries so no point
        // continuing
        return features;
      }
    }

    // System.out.println(startIndex + " " + endIndex);

    for (int i = startIndex; i <= endIndex; ++i) {
      features.add(locations.get(i));
    }

    return features;
  }
}
