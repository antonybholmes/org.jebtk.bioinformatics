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

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.jebtk.bioinformatics.genomic.Genome;

/**
 * Server for genome feature annotations.
 *
 * @author Antony Holmes
 *
 */
public class FeaturesServer {

  /**
   * The constant INSTANCE.
   */
  private static final FeaturesServer INSTANCE = new FeaturesServer();

  /**
   * The constant DEFAULT_FEATURES_FILE.
   */
  public static final File DEFAULT_FEATURES_FILE = new File("res/features.xml");

  /**
   * Gets the single instance of FeaturesServer.
   *
   * @return single instance of FeaturesServer
   */
  public static final FeaturesServer getInstance() {
    return INSTANCE;
  }

  /**
   * The feature map.
   */
  // genome, group, feature name
  private Map<Genome, Map<String, Map<String, AbstractFeaturesSearch>>> mFeatureMap = new HashMap<Genome, Map<String, Map<String, AbstractFeaturesSearch>>>();

  /**
   * Instantiates a new features server.
   */
  private FeaturesServer() {
    // do nothing
  }

  /**
   * Load xml.
   *
   * @param file the file
   */
  public final void loadXml(File file) {
    try {
      System.out.println("Parsing " + file.getAbsolutePath());

      SAXParserFactory factory = SAXParserFactory.newInstance();
      SAXParser saxParser = factory.newSAXParser();

      FeaturesXmlHandler handler = new FeaturesXmlHandler();

      saxParser.parse(file.getAbsolutePath(), handler);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Adds the common region binary search.
   *
   * @param genome      the genome
   * @param group       the group
   * @param name        the name
   * @param description the description
   * @param filePrefix  the file prefix
   * @throws FileNotFoundException the file not found exception
   */
  public final void addCommonRegionBinarySearch(Genome genome, String group, String name, String description,
      String filePrefix) throws FileNotFoundException {
    AbstractFeaturesSearch features = new FeaturesCommonRegionBinarySearch(name, description, filePrefix);

    add(genome, group, features);
  }

  /**
   * Adds the binary search.
   *
   * @param genome      the genome
   * @param group       the group
   * @param name        the name
   * @param description the description
   * @param file        the file
   * @throws FileNotFoundException the file not found exception
   */
  public final void addBinarySearch(Genome genome, String group, String name, String description, Path file)
      throws FileNotFoundException {
    AbstractFeaturesSearch features = new FeaturesBinarySearch(name, description, file);

    add(genome, group, features);
  }

  /**
   * Adds the basic search.
   *
   * @param genome      the genome
   * @param group       the group
   * @param name        the name
   * @param description the description
   * @param file        the file
   */
  public final void addBasicSearch(Genome genome, String group, String name, String description, Path file) {
    AbstractFeaturesSearch features = new FeaturesBasicSearch(name, description, file);

    add(genome, group, features);
  }

  /**
   * Adds the.
   *
   * @param genome   the genome
   * @param group    the group
   * @param features the features
   */
  public final void add(Genome genome, String group, AbstractFeaturesSearch features) {
    if (!mFeatureMap.containsKey(genome)) {
      mFeatureMap.put(genome, new HashMap<String, Map<String, AbstractFeaturesSearch>>());
    }

    if (!mFeatureMap.get(genome).containsKey(group)) {
      mFeatureMap.get(genome).put(group, new HashMap<String, AbstractFeaturesSearch>());
    }

    System.err.println("adding feature " + genome + " " + group + " " + features.getName());

    mFeatureMap.get(genome).get(group).put(features.getName(), features);
  }

  /**
   * Gets the.
   *
   * @param genome the genome
   * @param group  the group
   * @param name   the name
   * @return the abstract features search
   */
  public final AbstractFeaturesSearch get(Genome genome, String group, String name) {
    System.err.println("feature " + genome + " " + group + " " + name);

    if (!mFeatureMap.containsKey(genome)) {
      return null;
    }

    if (!mFeatureMap.get(genome).containsKey(group)) {
      return null;
    }

    return mFeatureMap.get(genome).get(group).get(name);
  }

  /**
   * Returns the collection of features associated with a genome.
   *
   * @param genome the genome
   * @param group  the group
   * @return the collection
   */
  public final Collection<AbstractFeaturesSearch> get(Genome genome, String group) {
    if (!mFeatureMap.containsKey(genome)) {
      return null;
    }

    if (!mFeatureMap.get(genome).containsKey(group)) {
      return null;
    }

    return mFeatureMap.get(genome).get(group).values();
  }

  /**
   * Uncaches any loaded feature sets.
   */
  public final void freeMemory() {
    for (Genome genome : mFeatureMap.keySet()) {
      for (String group : mFeatureMap.get(genome).keySet()) {
        for (String name : mFeatureMap.get(genome).get(group).keySet()) {
          mFeatureMap.get(genome).get(group).get(name).freeMemory();
        }
      }
    }
  }
}