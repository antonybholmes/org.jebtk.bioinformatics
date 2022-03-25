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
package org.jebtk.bioinformatics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * Server for gene synonyms.
 *
 * @author Antony Holmes
 */
public class GeneSynonymsService {

  /**
   * The constant INSTANCE.
   */
  private static final GeneSynonymsService INSTANCE = new GeneSynonymsService();

  /**
   * The constant DEFAULT_GENES_FILE.
   */
  public static final File DEFAULT_GENES_FILE = new File("res/gene_synonyms.txt");

  /**
   * Gets the single instance of GeneSynonymsService.
   *
   * @return single instance of GeneSynonymsService
   */
  public static final GeneSynonymsService getInstance() {
    return INSTANCE;
  }

  /**
   * The map.
   */
  // genome, group, feature name
  private Map<String, Set<String>> map = new HashMap<String, Set<String>>();

  /**
   * Instantiates a new gene synonyms service.
   */
  private GeneSynonymsService() {
    try {
      load();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Load.
   *
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public final void load() throws IOException {
    BufferedReader reader = new BufferedReader(new FileReader(DEFAULT_GENES_FILE));

    String line;
    List<String> tokens;

    try {
      reader.readLine();

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.fastSplit(line, TextUtils.TAB_DELIMITER);

        if (!map.containsKey(tokens.get(0))) {
          map.put(tokens.get(0), new HashSet<String>());
        }

        map.get(tokens.get(0)).add(tokens.get(1));
      }
    } finally {
      reader.close();
    }
  }

  /**
   * Gets the synonyms.
   *
   * @param name the name
   * @return the synonyms
   */
  public Set<String> getSynonyms(String name) {
    if (!map.containsKey(name)) {
      return null;
    }

    return new HashSet<String>(map.get(name));
  }
}