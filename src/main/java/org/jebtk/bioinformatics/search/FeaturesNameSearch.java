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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.ChromosomeService;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * The class FeaturesNameSearch.
 */
public class FeaturesNameSearch implements Iterable<Feature> {

  /**
   * The features.
   */
  protected Map<String, Feature> mFeatures = new HashMap<String, Feature>();

  /**
   * The type.
   */
  protected String mType;

  /**
   * Instantiates a new features name search.
   *
   * @param type the type
   */
  public FeaturesNameSearch(String type) {
    mType = type;
  }

  /**
   * Gets the type.
   *
   * @return the type
   */
  public final String getType() {
    return mType;
  }

  /**
   * Cache features.
   *
   * @param file the file
   */
  public void cacheFeatures(Path file) {
    mFeatures.clear();

    short chromosomeCount = 0;
    int featureCount = 0;

    try {
      BufferedReader reader = FileUtils.newBufferedReader(file);

      String line;

      try {
        // skip header
        line = reader.readLine();

        while ((line = reader.readLine()) != null) {
          if (Io.isEmptyLine(line)) {
            continue;
          }

          List<String> row = TextUtils.fastSplitRemoveQuotes(line);

          Feature feature = new Feature(row.get(0), ChromosomeService.getInstance().guessChr(file, row.get(1)),
              Integer.parseInt(row.get(2)), Integer.parseInt(row.get(3)));

          mFeatures.put(feature.getName(), feature);

          ++featureCount;
        }
      } finally {
        System.out
            .println("Loaded " + chromosomeCount + " groups with " + featureCount + " features from " + file + ".");

        reader.close();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Gets the feature.
   *
   * @param name the name
   * @return the feature
   */
  public final Feature getFeature(String name) {
    if (mFeatures.containsKey(name)) {
      return mFeatures.get(name);
    }

    return null;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  public Iterator<Feature> iterator() {
    return mFeatures.values().iterator();
  }

  /**
   * Size.
   *
   * @return the int
   */
  public int size() {
    return mFeatures.size();
  }
}
