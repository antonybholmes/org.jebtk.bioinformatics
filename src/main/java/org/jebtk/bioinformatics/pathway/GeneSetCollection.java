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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;

import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.model.ListModel;
import org.jebtk.core.text.TextUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Represents a gene along with a description of the gene,.
 *
 * @author Antony Holmes
 */
public class GeneSetCollection extends ListModel<GeneSet> implements Comparable<GeneSetCollection> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The member name.
   */
  private String mName;

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(GeneSetCollection.class);

  /**
   * Instantiates a new gene set collection.
   *
   * @param name the name
   */
  public GeneSetCollection(String name) {
    mName = name;
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public String getName() {
    return mName;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.util.AbstractList#hashCode()
   */
  @Override
  public int hashCode() {
    return mName.hashCode();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(GeneSetCollection o) {
    return mName.compareTo(o.mName);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.util.AbstractList#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof GeneSetCollection) {
      return compareTo((GeneSetCollection) o) == 0;
    } else {
      return false;
    }
  }

  /**
   * Gets the gene set name.
   *
   * @param file the file
   * @return the gene set name
   */
  public static String getGeneSetName(Path file) {
    return PathUtils.getName(file).toLowerCase().replaceFirst("\\..+", "");
  }

  /**
   * Parses the.
   *
   * @param file the file
   * @return the gene set collection
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static GeneSetCollection parse(Path file) throws IOException {
    LOG.info("Parsing gene set collection in {}...", file);

    return parseGMT(file);
  }

  /**
   * Parses the.
   *
   * @param is   the is
   * @param name the name
   * @return the gene set collection
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static GeneSetCollection parseGMT(Path file) throws IOException {
    String name = getGeneSetName(file);

    GeneSetCollection collection = new GeneSetCollection(name);

    BufferedReader reader = FileUtils.newBufferedReader(file);

    String line;
    List<String> tokens;

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.tabSplit(line);

        GeneSet geneSet = new GeneSet(tokens.get(0), name, CollectionUtils.subList(tokens, 2));

        collection.add(geneSet);
      }
    } finally {
      reader.close();
    }

    return collection;
  }

  /**
   * Parses the.
   *
   * @param file       the file
   * @param collection the collection
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void parse(Path file, Set<GeneSet> collection) throws IOException {
    LOG.info("Parsing gene set collection in {}...", file);

    String name = getGeneSetName(file);

    InputStream is = FileUtils.newBufferedInputStream(file);

    try {
      parse(is, name, collection);
    } finally {
      is.close();
    }
  }

  /**
   * Parses the.
   *
   * @param is         the is
   * @param name       the name
   * @param collection the collection
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void parse(InputStream is, String name, Set<GeneSet> collection) throws IOException {
    BufferedReader reader = new BufferedReader(new InputStreamReader(is));

    String line;
    List<String> tokens;

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.tabSplit(line);

        GeneSet geneSet = new GeneSet(tokens.get(0), name, CollectionUtils.subList(tokens, 2));

        collection.add(geneSet);
      }
    } finally {
      reader.close();
    }
  }
}
