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

import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultHashMap;
import org.jebtk.core.collections.EntryCreator;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.objectdb.RadixObjectDb;
import org.jebtk.core.objectdb.RadixObjectNode;
import org.jebtk.core.tree.TreeNode;

/**
 * The class MotifsFile.
 */
public class MotifsFs extends MotifDataSource implements Comparable<MotifsFs> {

  /**
   * The member dir.
   */
  protected Path mDir;

  /**
   * Instantiates a new motifs file.
   *
   * @param dir the dir
   */
  public MotifsFs(Path dir) {
    mDir = dir;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    return mDir.hashCode();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof MotifsFs) {
      return compareTo((MotifsFs) o) == 0;
    } else {
      return false;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(MotifsFs fs) {
    return mDir.compareTo(fs.mDir);
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

    // process dirs first

    List<Path> dirs = FileUtils.lsdir(root);

    // Only create a folder if it has something in it

    if (terms != null && terms.size() > 0) {

      // For each dir, create a radix db for searching
      Map<Path, RadixObjectDb<Motif>> motifMap = DefaultHashMap.create(new EntryCreator<RadixObjectDb<Motif>>() {
        @Override
        public RadixObjectDb<Motif> newEntry() {
          return new RadixObjectDb<Motif>();
        }
      });

      Map<Path, TreeNode<Motif>> nodeMap = new HashMap<Path, TreeNode<Motif>>();

      for (Path dir : dirs) {
        nodeMap.put(dir, new TreeNode<Motif>(PathUtils.getName(dir)));

        if (FileUtils.ls(dir).size() > 0) {
          List<Path> files = FileUtils.ls(dir);

          for (Path file : files) {
            if (PathUtils.getName(file).endsWith("motif.gz")) {
              // Extract all the nodes from the file and add them
              for (Motif motif : Motif.parseMotifs(file, PathUtils.getName(file.getParent()))) {

                if (caseSensitive) {
                  motifMap.get(dir).addObject(motif.getName(), motif);
                  motifMap.get(dir).addObject(motif.getId(), motif);
                  motifMap.get(dir).addObject(motif.getGene(), motif);
                } else {
                  motifMap.get(dir).addObject(motif.getName().toLowerCase(), motif);
                  motifMap.get(dir).addObject(motif.getId().toLowerCase(), motif);
                  motifMap.get(dir).addObject(motif.getGene().toLowerCase(), motif);
                }
              }
            }
          }
        }
      }

      for (String term : terms) {
        String s = caseSensitive ? term : term.toLowerCase();

        for (Path dir : dirs) {
          if (motifMap.containsKey(dir)) {

            RadixObjectNode<Motif> db = motifMap.get(dir).getChild(s);

            if (db != null) {
              for (Motif motif : db) {
                nodeMap.get(dir).addChild(new TreeNode<Motif>(motif.getName() + " (" + motif.getId() + ")", motif));
              }
            }
          }
        }
      }

      for (Path dir : dirs) {
        TreeNode<Motif> node = nodeMap.get(dir);

        // Only add dir (database) nodes if they have some children
        if (node.getChildCount() > 0) {
          rootNode.addChild(node);
        }
      }
    } else {
      // Add all nodes if there are no search terms.

      for (Path dir : dirs) {
        if (FileUtils.ls(dir).size() > 0) {
          TreeNode<Motif> node = new TreeNode<Motif>(PathUtils.getName(dir));

          List<Path> files = FileUtils.ls(dir);

          for (Path file : files) {
            if (PathUtils.getName(file).endsWith("motif.gz")) {
              // Extract all the nodes from the file and add them
              for (Motif motif : CollectionUtils.sort(Motif.parseMotifs(file, PathUtils.getName(file.getParent())))) {
                node.addChild(new TreeNode<Motif>(motif.getName() + " (" + motif.getId() + ")", motif));
              }
            }
          }

          if (node.getChildCount() > 0) {
            rootNode.addChild(node);
          }
        }
      }
    }
  }

  protected static void filter(Motifs motifs, TreeNode<Motif> rootNode, List<String> terms, boolean inList,
      boolean exactMatch, boolean caseSensitive) throws Exception {

    TreeNode<Motif> node = new TreeNode<Motif>(motifs.getName());

    int count = 0;

    for (Motif motif : motifs) {
      if (terms.size() == 0) {
        node.addChild(new TreeNode<Motif>(motif.getName() + " (" + motif.getId() + ")", motif));
        ++count;
        continue;
      }

      boolean found = false;

      for (String term : terms) {
        if (caseSensitive) {
          if (exactMatch) {
            if (motif.getName().equals(term) || motif.getId().equals(term) || motif.getGene().equals(term)) {
              found = true;
              break;
            }
          } else {
            if (motif.getName().contains(term) || motif.getId().contains(term) || motif.getGene().contains(term)) {
              found = true;
              break;
            }
          }
        } else {
          String lcs = term.toLowerCase();

          if (exactMatch) {
            if (motif.getName().toLowerCase().equals(lcs) || motif.getId().toLowerCase().equals(lcs)
                || motif.getGene().toLowerCase().equals(lcs)) {
              found = true;
              break;
            }
          } else {
            if (motif.getName().toLowerCase().contains(lcs) || motif.getId().toLowerCase().contains(lcs)
                || motif.getGene().toLowerCase().contains(lcs)) {
              found = true;
              break;
            }
          }
        }
      }

      if ((found && inList) || (!found && !inList)) {
        node.addChild(new TreeNode<Motif>(motif.getName() + " (" + motif.getId() + ")", motif));
        ++count;
      }
    }

    if (count > 0) {
      rootNode.addChild(node);
    }
  }
}
