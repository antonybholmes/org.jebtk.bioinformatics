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
import java.io.BufferedWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.text.TextUtils;

/**
 * The class FeaturesCommonRegionBinarySearchWriter.
 */
public class FeaturesCommonRegionBinarySearchWriter extends FeaturesBasicSearch {

  /**
   * Instantiates a new features common region binary search writer.
   *
   * @param name the name
   * @param file the file
   */
  public FeaturesCommonRegionBinarySearchWriter(String name, Path file) {
    super(name, null, file);

    cacheFeatures();
  }

  /**
   * Write.
   */
  public void write() {
    Path dir = mFile.getParent();

    Path locationsPath = null;

    locationsPath = dir.resolve(mName + "_bins_locations.txt");

    Path binsPath = null;

    binsPath = dir.resolve(mName + "_bins.txt");

    System.out.println("writing " + locationsPath);
    System.out.println("writing " + binsPath);

    String header = "";

    try {
      BufferedReader reader = FileUtils.newBufferedReader(mFile);

      try {
        header = reader.readLine();
      } finally {
        reader.close();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

    header += TextUtils.TAB_DELIMITER + "bins";

    try {
      BufferedWriter locationsWriter = FileUtils.newBufferedWriter(dir.resolve(locationsPath));

      BufferedWriter featuresWriter = FileUtils.newBufferedWriter(dir.resolve(binsPath));

      try {
        featuresWriter.write(header);
        featuresWriter.newLine();

        locationsWriter.write("chromosome");
        locationsWriter.write(TextUtils.TAB_DELIMITER);
        locationsWriter.write("location");
        locationsWriter.newLine();

        for (short chromosome = 1; chromosome < allLocations.size(); ++chromosome) {
          if (allLocations.get(chromosome) == null) {
            continue;
          }

          List<Integer> locations = new ArrayList<Integer>();

          Set<Integer> inUse = new HashSet<Integer>();

          for (Feature feature : allLocations.get(chromosome)) {
            if (!inUse.contains(feature.getStart())) {
              locations.add(feature.getStart());
              inUse.add(feature.getStart());
            }

            if (!inUse.contains(feature.getEnd())) {
              locations.add(feature.getEnd());
              inUse.add(feature.getEnd());
            }
          }

          // sort them

          System.out.println("chr" + chromosome);

          Collections.sort(locations);

          for (int location : locations) {
            locationsWriter.write(Integer.toString(chromosome));
            locationsWriter.write(TextUtils.TAB_DELIMITER);
            locationsWriter.write(Integer.toString(location));
            locationsWriter.newLine();
          }

          // see which features overlap the start of a location
          for (Feature feature : allLocations.get(chromosome)) {
            List<Integer> overlap = new ArrayList<Integer>();

            for (int i = 0; i < locations.size(); ++i) {
              int location = locations.get(i);

              if (feature.getStart() == location || (feature.getStart() < location && feature.getEnd() > location)) {
                // every feature will be added to the bin where it starts
                // also add features that overlap into the bin

                // Each bin is searched from the start to the beginning of the
                // next bin (not including the end position, which is the start
                // of
                // the next bin) so only features that start where this bin does
                // or overlap it may be added otherwise the feature
                // is allocated to another bin
                overlap.add(i);
              }
            }

            featuresWriter.write(feature.getName());
            featuresWriter.write(TextUtils.TAB_DELIMITER);
            featuresWriter.write(feature.getChr().toString());
            featuresWriter.write(TextUtils.TAB_DELIMITER);
            featuresWriter.write(Integer.toString(feature.getStart()));
            featuresWriter.write(TextUtils.TAB_DELIMITER);
            featuresWriter.write(Integer.toString(feature.getEnd()));
            featuresWriter.write(TextUtils.TAB_DELIMITER);
            featuresWriter.write(TextUtils.join(overlap, TextUtils.COMMA_DELIMITER));
            featuresWriter.newLine();
          }
        }
      } finally {
        locationsWriter.close();
        featuresWriter.close();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
