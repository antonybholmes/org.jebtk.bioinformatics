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

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.jebtk.core.Resources;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.tree.TreeNode;
import org.xml.sax.SAXException;

/**
 * The class MotifsFile.
 */
public class MotifsXmlFs extends MotifsFs {

  /**
   * Instantiates a new motifs file.
   *
   * @param dir the dir
   */
  public MotifsXmlFs(Path dir) {
    super(dir);
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.motifs.MotifsFs#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof MotifsXmlFs) {
      return compareTo((MotifsXmlFs) o) == 0;
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
    if (!FileUtils.exists(root)) {
      return;
    }

    List<Path> files = FileUtils.ls(root, false, true);

    for (Path file : files) {
      if (!PathUtils.getName(file).endsWith("xml.gz")) {
        continue;
      }

      Motifs motifs = parseMotifXmlGz(file);

      filter(motifs, rootNode, terms, inList, exactMatch, caseSensitive);
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
  public static Motifs parseMotifXml(Path file) throws IOException, ParserConfigurationException, SAXException {
    InputStream stream = FileUtils.newBufferedInputStream(file);

    Motifs motifs = null;

    try {
      motifs = parseMotifXml(stream);
    } finally {
      stream.close();
    }

    return motifs;
  }

  /**
   * Parses the motif xml gz.
   *
   * @param file the file
   * @return the motifs
   * @throws IOException                  Signals that an I/O exception has
   *                                      occurred.
   * @throws ParserConfigurationException the parser configuration exception
   * @throws SAXException                 the SAX exception
   */
  public static Motifs parseMotifXmlGz(Path file) throws IOException, ParserConfigurationException, SAXException {
    InputStream stream = Resources.getGzipInputStream(file);

    Motifs motifs = null;

    try {
      motifs = parseMotifXml(stream);
    } finally {
      stream.close();
    }

    return motifs;
  }

  /**
   * Parses the motif xml.
   *
   * @param is the is
   * @return the motifs
   * @throws ParserConfigurationException the parser configuration exception
   * @throws SAXException                 the SAX exception
   * @throws IOException                  Signals that an I/O exception has
   *                                      occurred.
   */
  public static Motifs parseMotifXml(InputStream is) throws ParserConfigurationException, SAXException, IOException {
    if (is == null) {
      return null;
    }

    SAXParserFactory factory = SAXParserFactory.newInstance();
    SAXParser saxParser = factory.newSAXParser();

    MotifXmlHandler handler = new MotifXmlHandler();

    saxParser.parse(is, handler);

    return handler.getMotifs();
  }
}
