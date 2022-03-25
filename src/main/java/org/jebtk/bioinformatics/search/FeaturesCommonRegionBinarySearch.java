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
import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArrayList;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.ChromosomeService;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.settings.SettingsService;
import org.jebtk.core.text.TextUtils;

/**
 * The class FeaturesCommonRegionBinarySearch.
 */
public class FeaturesCommonRegionBinarySearch extends AbstractFeaturesSearch {

  /**
   * The constant BIN_FEATURE_FILE_ENDING.
   */
  public static final String BIN_FEATURE_FILE_ENDING = "_bins.txt";

  /**
   * The constant BIN_LOCATION_FILE_ENDING.
   */
  public static final String BIN_LOCATION_FILE_ENDING = "_bins_locations.txt";

  /**
   * The features at location by index.
   */
  private List<List<FeatureBin>> mFeaturesAtLocationByIndex = new CopyOnWriteArrayList<List<FeatureBin>>();

  /**
   * The bins file.
   */
  private Path mBinsFile;

  /**
   * The feature file.
   */
  private Path mFeatureFile;

  /**
   * The size.
   */
  private int size = 0;

  /**
   * Instantiates a new features common region binary search.
   *
   * @param name        the name
   * @param description the description
   * @param filePrefix  the file prefix
   * @throws FileNotFoundException the file not found exception
   */
  public FeaturesCommonRegionBinarySearch(String name, String description, String filePrefix)
      throws FileNotFoundException {
    super(name, description);

    mBinsFile = PathUtils.getPath(filePrefix + BIN_LOCATION_FILE_ENDING);
    mFeatureFile = PathUtils.getPath(filePrefix + BIN_FEATURE_FILE_ENDING);

    System.err.println(
        "Feature prefix " + filePrefix + " " + FileUtils.exists(mBinsFile) + " " + FileUtils.exists(mFeatureFile));

  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return PathUtils.toString(mFeatureFile);
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * cacheFeatures()
   */
  public final void cacheFeatures() {
    mFeaturesAtLocationByIndex.clear();

    Chromosome chromosome;
    int startLocation;

    // first the bins file

    // since we have 25 chromosomes, create place holders for them
    // we need one more position than necessary because index 0 is
    // never used since chromosomes are one based indexes.

    for (int i = 0; i < 26; ++i) {
      mFeaturesAtLocationByIndex.add(new ArrayList<FeatureBin>());
    }

    try {
      BufferedReader reader = FileUtils.newBufferedReader(mBinsFile);

      String line;

      try {
        // skip header
        reader.readLine();

        while ((line = reader.readLine()) != null) {
          if (Io.isEmptyLine(line)) {
            continue;
          }

          List<String> row = TextUtils.fastSplit(line, TextUtils.TAB_DELIMITER);

          chromosome = ChromosomeService.getInstance().guessChr(mFeatureFile, row.get(0));

          startLocation = Integer.parseInt(row.get(1));

          // System.out.println("adding feature " + type + " " + chromosome + "
          // " +
          // startLocation + " " + endLocation);

          // start

          mFeaturesAtLocationByIndex.get(chromosome.getId()).add(new FeatureBin(startLocation));

        }
      } finally {
        reader.close();
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

    // now the features

    try {
      BufferedReader reader = FileUtils.newBufferedReader(mFeatureFile);

      String line;

      try {
        // skip header
        reader.readLine();

        while ((line = reader.readLine()) != null) {
          if (Io.isEmptyLine(line)) {
            continue;
          }

          List<String> row = TextUtils.tabSplit(line);

          Feature feature = new Feature(row.get(0), ChromosomeService.getInstance().guessChr(mFeatureFile, row.get(1)),
              Integer.parseInt(row.get(2)), Integer.parseInt(row.get(3)));

          List<String> bins = TextUtils.fastSplit(row.get(4), TextUtils.COMMA_DELIMITER);

          for (String bin : bins) {
            int b = Integer.parseInt(bin);

            // System.out.println(feature.getChromosome() + " " + b + " " +
            // featuresAtLocationByIndex.get(feature.getChromosome()).size() + "
            // " +
            // featuresAtLocationByIndex.get(feature.getChromosome()).get(featuresAtLocationByIndex.get(feature.getChromosome()).size()
            // - 1).getStart());

            // add the feature to each of the bins it belongs to
            mFeaturesAtLocationByIndex.get(feature.getChr().getId()).get(b).add(feature);

            ++size;
          }

          // System.out.println(feature.getName());

          mFeatureByName.put(feature.getName(), feature);
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
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome, int, int)
   */
  @Override
  public List<Feature> getFeatures(Chromosome chromosome, int startLocation, int endLocation) {
    if (mFeaturesAtLocationByIndex.size() == 0) {
      cacheFeatures();
    }

    int startIndex;
    int endIndex;

    List<Feature> returnFeatures = new ArrayList<Feature>();

    List<FeatureBin> features = mFeaturesAtLocationByIndex.get(chromosome.getId());

    if (features.size() == 0) {
      return returnFeatures;
    }

    // System.out.println(binsFile + ":" + featuresAtLocationByIndex.size() +
    // ":" +
    // chromosome);

    int last = features.size() - 1;

    if (startLocation > features.get(last).getStart() || endLocation < features.get(0).getStart()) {
      // the range is clearly not within the feature set so don't even bother to
      // look
      return returnFeatures;
    }

    if (startLocation <= features.get(0).getStart()) {
      startIndex = 0;
    } else {
      startIndex = binarySearchInner(startLocation, chromosome, features);
    }

    if (endLocation >= features.get(last).getStart()) {
      endIndex = last;
    } else {
      endIndex = binarySearchOuter(endLocation, chromosome, features);
    }

    if (startIndex > endIndex) {
      // this should never happen
      return returnFeatures;
    }

    // System.out.println(startIndex + " " + endIndex + " " + chromosome + " "
    // +startLocation + " " + features.get(0).getStart() + " " +
    // features.get(last).getStart() + " " + features.size());

    for (int i = startIndex; i <= endIndex; ++i) {
      // add all the features
      for (Feature feature : features.get(i)) {
        returnFeatures.add(feature);
      }
    }

    return returnFeatures;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome, int)
   */
  @Override
  public List<Feature> getFeatures(Chromosome chromosome, int location) {
    if (mFeaturesAtLocationByIndex.size() == 0) {
      cacheFeatures();
    }

    List<Feature> returnFeatures = new ArrayList<Feature>();

    List<FeatureBin> features = mFeaturesAtLocationByIndex.get(chromosome.getId());

    if (features.size() == 0) {
      return returnFeatures;
    }

    // System.out.println(binsFile + ":" + featuresAtLocationByIndex.size() +
    // ":" +
    // chromosome);

    int index = binarySearch(location, chromosome, features);

    if (index == -1) {
      return returnFeatures;
    }

    for (Feature feature : features.get(index)) {
      returnFeatures.add(feature);
    }

    return returnFeatures;
  }

  /**
   * Binary search inner.
   *
   * @param value      the value
   * @param chromosome the chromosome
   * @param features   the features
   * @return the int
   */
  private static final int binarySearchInner(int value, Chromosome chromosome, List<FeatureBin> features) {
    int location = binarySearch(value, chromosome, features);

    // System.out.println("loc" + location);

    if (features.size() > 0) {
      return location;
    }

    // we are in a gap between features so return the next one along
    return location + 1;
  }

  /**
   * Binary search outer.
   *
   * @param value      the value
   * @param chromosome the chromosome
   * @param features   the features
   * @return the int
   */
  private static final int binarySearchOuter(int value, Chromosome chromosome, List<FeatureBin> features) {
    int location = binarySearch(value, chromosome, features);

    if (features.size() > 0) {
      return location;
    }

    // we are in a gap between features so return the next one along
    return location - 1;
  }

  /**
   * Binary search.
   *
   * @param value      the value
   * @param chromosome the chromosome
   * @param features   the features
   * @return the int
   */
  private static final int binarySearch(int value, Chromosome chromosome, List<FeatureBin> features) {
    int min = 0;
    int max = features.size() - 1;

    int mid = -1;

    int midMax = max;

    while (min <= max) {
      mid = (min + max) / 2;

      // System.out.println("s:" + value + " " + mid + " " + midMax + " min:" +
      // min +
      // " max:" + max + " " + locations.get(mid).getStart() + " " +
      // locations.get(mid
      // + 1).getStart());

      if (mid == midMax) {
        // special case where value falls exactly in the last bin so we
        // check if it occurs in the last bin
        if (features.get(mid).getStart() <= value) {
          return mid;
        }

        // does not occur in the last bin which will trigger
        // -1 to be returned. In a correctly formed expression and
        // feature file, this should never happen, but we allow for
        // it as a check for quick debugging purposes.
        break;
      }

      if (features.get(mid).getStart() <= value && features.get(mid + 1).getStart() > value) {
        return mid;
      } else if (value >= features.get(mid + 1).getStart()) {
        min = mid + 1;
      } else {
        max = mid - 1;
      }
    }

    // we should not return -1 unless the coordinate
    // is completely out of the chromosome bounds, but this
    // is checked for in a previous step. We should only be
    // searching within bounds.
    return -1;
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * freeMemory( )
   */
  public void freeMemory() {
    super.freeMemory();

    mFeaturesAtLocationByIndex.clear();
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.search.AbstractFeaturesSearch#
   * getFeatures (edu.columbia.rdf.lib.bioinformatics.genome.Chromosome)
   */
  public List<Feature> getFeatures(Chromosome chromosome) {
    if (mFeaturesAtLocationByIndex.size() == 0) {
      cacheFeatures();
    }

    List<Feature> features = new ArrayList<Feature>();

    Set<String> used = new HashSet<String>();

    for (FeatureBin featureBin : mFeaturesAtLocationByIndex.get(chromosome.getId())) {
      for (Feature feature : featureBin) {
        if (used.contains(feature.toString())) {
          continue;
        }

        used.add(feature.toString());
        features.add(feature);
      }
    }

    return features;
  }

  /**
   * Returns the feature file associated with the database genomic feature name.
   *
   * @param group the group
   * @param name  the name
   * @return the features file
   */
  public static final File getFeaturesFile(String group, String name) {
    return new File(SettingsService.getInstance().getString(group + "." + name) + BIN_FEATURE_FILE_ENDING);
  }

  /**
   * Returns the location file associated with the genomic feature name.
   *
   * @param group the group
   * @param name  the name
   * @return the locations file
   */
  public static final File getLocationsFile(String group, String name) {
    return new File(SettingsService.getInstance().getString(group + "." + name) + BIN_LOCATION_FILE_ENDING);
  }

  /**
   * Gets the feature file.
   *
   * @return the feature file
   */
  public Path getFeatureFile() {
    return mFeatureFile;
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
