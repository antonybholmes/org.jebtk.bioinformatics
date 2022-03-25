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
package org.jebtk.bioinformatics.motifs;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.jebtk.bioinformatics.BaseCounts;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.TextUtils;
import org.jebtk.core.tree.TreeNode;
import org.xml.sax.SAXException;

/**
 * The class MotifsFile.
 */
public class MotifsPwtFs extends MotifsFs {

  /**
   * Instantiates a new motifs file.
   *
   * @param dir the dir
   */
  public MotifsPwtFs(Path dir) {
    super(dir);
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.motifs.MotifsFs#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof MotifsPwtFs) {
      return compareTo((MotifsPwtFs) o) == 0;
    } else {
      return false;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.motifs.MotifsDB#createTree(org.abh.lib.
   * tree.TreeRootNode, java.lang.String)
   */
  @Override
  public void createTree(TreeNode<Motif> root, List<String> terms, boolean inList, boolean exactMatch,
      boolean caseSensitive) throws Exception {
    // TreeRootNode<Motif> root = new TreeRootNode<Motif>();

    createTreeDir(mDir, root, terms, inList, exactMatch, caseSensitive);

    // return root;
  }

  /**
   * Creates the tree dir.
   *
   * @param root          the root
   * @param rootNode      the root node
   * @param terms         the terms
   * @param inList        the in list
   * @param exactMatch    the exact match
   * @param caseSensitive the case sensitive
   * @return the int
   * @throws Exception the exception
   */
  private static void createTreeDir(Path root, TreeNode<Motif> rootNode, List<String> terms, boolean inList,
      boolean exactMatch, boolean caseSensitive) throws Exception {
    System.err.println("root " + root);
    if (!FileUtils.exists(root)) {
      return;
    }

    List<Path> files = FileUtils.ls(root, false, true);

    for (Path file : files) {
      if (PathUtils.getName(file).endsWith("pwt.gz")) {
        Motifs motifs = parseMotifPwt(file);

        filter(motifs, rootNode, terms, inList, exactMatch, caseSensitive);
      }
    }
  }

  /**
   * Parses the motif xml.
   *
   * @param file the file
   * @return the motifs
   * @throws IOException                  Signals that an I/O exception has
   *                                      occurred.
   * @throws ParserConfigurationException the parser configuration exception
   * @throws SAXException                 the SAX exception
   */
  public static Motifs parseMotifPwt(Path file) throws IOException, ParserConfigurationException, SAXException {
    BufferedReader is = FileUtils.newBufferedReader(file);

    List<Motif> motifs = new ArrayList<Motif>();

    String db = null;

    try {
      // Skip 'db:' prefix
      db = is.readLine().substring(3);

      // Skip header
      is.readLine();

      String line;

      while ((line = is.readLine()) != null) {
        List<String> tokens = TextUtils.tabSplit(line);

        String id = tokens.get(0);
        String name = tokens.get(1);
        int l = Integer.parseInt(tokens.get(2));
        List<String> bases = TextUtils.scSplit(tokens.get(3));

        List<BaseCounts> counts = new ArrayList<BaseCounts>(l);

        for (String base : bases) {
          List<String> values = TextUtils.commaSplit(base);

          double a = Double.parseDouble(values.get(0));
          double c = Double.parseDouble(values.get(1));
          double g = Double.parseDouble(values.get(2));
          double t = Double.parseDouble(values.get(3));

          counts.add(new BaseCounts(a, c, g, t, true));
        }

        motifs.add(new Motif(id, name, name, db, counts));
      }

    } finally {
      is.close();
    }

    return new Motifs(db, motifs);
  }
}
