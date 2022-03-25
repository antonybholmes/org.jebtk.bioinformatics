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

import java.util.List;

import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.tree.TreeNode;

/**
 * The class MotifsDBService.
 */
public class MotifsDataSourceService extends MotifDataSource {

  /**
   * The Class MotifsDBServiceLoader.
   */
  private static class MotifsDBServiceLoader {

    /** The Constant INSTANCE. */
    private static final MotifsDataSourceService INSTANCE = new MotifsDataSourceService();
  }

  /**
   * Gets the single instance of SettingsService.
   *
   * @return single instance of SettingsService
   */
  public static MotifsDataSourceService getInstance() {
    return MotifsDBServiceLoader.INSTANCE;
  }

  /**
   * The member dbs. We only load each db type once
   */
  private List<MotifDataSource> mDbs = new UniqueArrayList<MotifDataSource>();

  /**
   * Instantiates a new motifs db service.
   */
  private MotifsDataSourceService() {
    // do nothing
  }

  /**
   * Adds the back end.
   *
   * @param motifsDb the motifs db
   */
  public void addDataSource(MotifDataSource motifsDb) {
    mDbs.add(motifsDb);
  }

  /**
   * Gets the DB count.
   *
   * @return the DB count
   */
  public int getDBCount() {
    return mDbs.size();
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

    for (MotifDataSource db : mDbs) {
      db.createTree(root, terms, inList, exactMatch, caseSensitive);
    }

    // Sort the folders
    root.sortChildren();
  }
}
