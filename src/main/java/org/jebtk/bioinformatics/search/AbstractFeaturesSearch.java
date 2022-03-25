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

import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.bioinformatics.Bounds;
import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;

/**
 * The class AbstractFeaturesSearch.
 */
public abstract class AbstractFeaturesSearch {

  /**
   * The feature by name.
   */
  protected Map<String, Feature> mFeatureByName = new HashMap<String, Feature>();

  /**
   * The member name.
   */
  protected String mName;

  /**
   * The member description.
   */
  private String mDescription;

  /**
   * Instantiates a new abstract features search.
   *
   * @param name        the name
   * @param description the description
   */
  public AbstractFeaturesSearch(String name, String description) {
    mName = name;
    mDescription = description;
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public final String getName() {
    return mName;
  }

  /**
   * Should return description or common name of feature set if it is not
   * described in the name.
   *
   * @return the description
   */
  public final String getDescription() {
    return mDescription;
  }

  /*
   * public void cacheFeaturesBin(File file) { allLocations.clear();
   * 
   * // number of items
   * 
   * int stringSize;
   * 
   * short chromosome;
   * 
   * int startLocation; int endLocation;
   * 
   * String name;
   * 
   * short chromosomeCount = 0; int featureCount = 0;
   * 
   * try { DataInputStream dis = new DataInputStream(new BufferedInputStream(new
   * FileInputStream(file)));
   * 
   * try { //size dis.readInt();
   * 
   * stringSize = dis.read();
   * 
   * //System.out.println("string size:" + stringSize);
   * 
   * while (true) { byte[] s = new byte[stringSize];
   * 
   * dis.read(s, 0, stringSize);
   * 
   * name = new String(s); name = name.trim();
   * 
   * chromosome = dis.read();
   * 
   * startLocation = dis.readInt();
   * 
   * endLocation = dis.readInt();
   * 
   * if (!allLocations.containsKey(chromosome)) { allLocations.put(chromosome, new
   * ArrayList<Feature>());
   * 
   * ++chromosomeCount; }
   * 
   * Feature feature = new Feature();
   * 
   * feature.getName() = name; feature.group = chromosome; feature.getStart() =
   * startLocation; feature.getEnd() = endLocation; //feature.type = type;
   * 
   * //System.out.println(name + " " + startLocation);
   * 
   * 
   * allLocations.get(chromosome).add(feature);
   * 
   * this.featureByName.put(feature.getName(), feature);
   * 
   * ++featureCount; } } finally { dis.close(); } } catch (Exception e) {
   * e.printStackTrace(); }
   * 
   * System.out.println("Loaded " + chromosomeCount + " groups with " +
   * featureCount + " features from " + file + "."); }
   */

  /**
   * Gets the names.
   *
   * @param chromosome    the chromosome
   * @param startLocation the start location
   * @param endLocation   the end location
   * @return the names
   * @throws ParseException the parse exception
   */
  public final Set<String> getNames(Chromosome chromosome, int startLocation, int endLocation) throws ParseException {
    List<Feature> features = getFeatures(chromosome, startLocation, endLocation);

    Set<String> names = new HashSet<String>();

    for (Feature feature : features) {
      names.add(feature.getName());
    }

    return names;
  }

  /**
   * Gets the sorted name list.
   *
   * @param chromosome    the chromosome
   * @param startLocation the start location
   * @param endLocation   the end location
   * @return the sorted name list
   * @throws ParseException the parse exception
   */
  public final List<String> getSortedNameList(Chromosome chromosome, int startLocation, int endLocation) {
    // go through the table

    List<String> names = getNameListInOrder(chromosome, startLocation, endLocation);

    Collections.sort(names);

    return names;
  }

  /**
   * Returns a unique list of feature names as they appear in the order they are
   * found (duplicates are removed).
   *
   * @param chromosome the chromosome
   * @param start      the start
   * @param end        the end
   * @return the name list in order
   * @throws ParseException the parse exception
   */
  public final List<String> getNameListInOrder(Chromosome chromosome, int start, int end) {
    // go through the table

    List<Feature> features = getFeatures(chromosome, start, end);

    List<String> names = new ArrayList<String>();

    Set<String> used = new HashSet<String>();

    for (Feature feature : features) {
      if (used.contains(feature.getName())) {
        continue;
      }

      names.add(feature.getName());
      used.add(feature.getName());
    }

    return names;
  }

  /**
   * Gets the feature.
   *
   * @param name the name
   * @return the feature
   */
  public final Feature getFeature(String name) {
    if (mFeatureByName.size() == 0) {
      cacheFeatures();
    }

    if (mFeatureByName.containsKey(name)) {
      return mFeatureByName.get(name);
    }

    return null;
  }

  /**
   * Gets the features.
   *
   * @param s the s
   * @return the features
   */
  public List<Feature> getFeatures(String s) {
    List<Feature> features = new ArrayList<Feature>();

    for (String name : this.mFeatureByName.keySet()) {
      if (name.toLowerCase().indexOf(s) != -1) {
        features.add(this.mFeatureByName.get(name));
      }
    }

    return features;
  }

  /**
   * Frees cached items from memory.
   */
  public void freeMemory() {
    System.out.println("Unloading " + this.mName);

    mFeatureByName.clear();
  }

  /**
   * Cache features.
   */
  public abstract void cacheFeatures();

  /**
   * Gets the features.
   *
   * @param chromosome the chromosome
   * @return the features
   * @throws ParseException the parse exception
   */
  public abstract List<Feature> getFeatures(Chromosome chromosome);

  public List<Feature> getFeatures(GenomicRegion region) {
    return getFeatures(region.getChr(), region.getStart(), region.getEnd());
  }

  /**
   * Gets the features.
   *
   * @param chromosome    the chromosome
   * @param startLocation the start location
   * @param endLocation   the end location
   * @return the features
   * @throws ParseException the parse exception
   */
  public abstract List<Feature> getFeatures(Chromosome chromosome, int start, int end);

  /**
   * Returns the set of features overlapping a location.
   *
   * @param chromosome the chromosome
   * @param location   the location
   * @return the features
   * @throws ParseException the parse exception
   */
  public abstract List<Feature> getFeatures(Chromosome chromosome, int location);

  /**
   * Returns the min and max bound from a list of features. This is useful for
   * ordered feature lists where the starts may be ordered but the ends are not so
   * the last element of the list does not necessarily contain the max end
   * location.
   *
   * @param features the features
   * @return the feature bounds
   */
  public static final Bounds getFeatureBounds(List<Feature> features) {
    Bounds bounds = new Bounds();

    for (Feature feature : features) {
      if (feature.getStart() < bounds.min) {
        bounds.min = feature.getStart();
      }

      if (feature.getEnd() > bounds.max) {
        bounds.max = feature.getEnd();
      }
    }

    return bounds;
  }

  /**
   * Produces a condensed string representation of cytoband information.
   *
   * @param features the features
   * @return a condensed string representation of cytoband information or an empty
   *         string if the list is empty/invalid.
   */
  public static final String condenseFeatures(List<Feature> features) {
    if (features == null || features.size() == 0) {
      return "";
    }

    if (features.size() == 1 || features.get(0).equals(features.get(features.size() - 1))) {
      return features.get(0).getName();
    } else {
      return features.get(0).getName() + "--" + features.get(features.size() - 1).getName();
    }
  }

  /**
   * Produces a condensed string representation of cytoband information.
   *
   * @param features the features
   * @return a condensed string representation of cytoband information or an empty
   *         string if the list is empty/invalid.
   */
  public static final String condense(List<String> features) {
    if (features == null || features.size() == 0) {
      return "";
    }

    if (features.size() == 1 || features.get(0).equals(features.get(features.size() - 1))) {
      return features.get(0);
    } else {
      return features.get(0) + "--" + features.get(features.size() - 1);
    }
  }

  /**
   * Size.
   *
   * @return the int
   */
  public abstract int size();
}
