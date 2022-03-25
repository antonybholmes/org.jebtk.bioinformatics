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

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.jebtk.core.model.SetModel;

/**
 * The class GeneSetCollections.
 */
public class GeneSetCollections extends SetModel<GeneSetCollection> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * Parses the.
   *
   * @param files the files
   * @return the gene set collections
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static GeneSetCollections parse(Path[] files) throws IOException {
    return parse(Arrays.asList(files));
  }

  /**
   * Parses the.
   *
   * @param files the files
   * @return the gene set collections
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static GeneSetCollections parse(List<Path> files) throws IOException {
    GeneSetCollections model = new GeneSetCollections();

    for (Path file : files) {
      GeneSetCollection collection = GeneSetCollection.parse(file);

      model.add(collection);
    }

    return model;
  }

  /**
   * Parses the.
   *
   * @param files       the files
   * @param collections the collections
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void parse(List<Path> files, Set<GeneSet> collections) throws IOException {
    for (Path file : files) {
      GeneSetCollection collection = GeneSetCollection.parse(file);

      for (GeneSet geneset : collection) {
        collections.add(geneset);
      }
    }
  }

  /*
   * public static GeneSetCollections loadFromResources() throws IOException {
   * GeneSetCollections model = new GeneSetCollections();
   * 
   * for (String res : Resources.getInstance()) { if (!res.contains(".gmt")) {
   * continue; }
   * 
   * String name = res.toLowerCase().replaceFirst("^.+\\/", "")
   * .replaceFirst("\\.symbols.+", "");
   * 
   * GeneSetCollection collection = GeneSetCollection
   * .parseGMT(Resources.getResInputStream(res), name);
   * 
   * model.add(collection); }
   * 
   * return model; }
   */
}
