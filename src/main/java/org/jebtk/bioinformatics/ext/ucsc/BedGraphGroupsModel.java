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
package org.jebtk.bioinformatics.ext.ucsc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jebtk.core.event.ChangeEvent;
import org.jebtk.core.event.ChangeListener;
import org.jebtk.core.event.ChangeListeners;

/**
 * A model for shared bed graphs.
 * 
 * @author Antony Holmes
 *
 */
public class BedGraphGroupsModel extends ChangeListeners implements ChangeListener, Iterable<BedGraphGroupModel> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant DEFAULT_GROUP.
   */
  private static final String DEFAULT_GROUP = "Default Group";

  /**
   * The member group map.
   */
  private Map<String, Integer> mGroupMap = new HashMap<String, Integer>();

  /**
   * The member groups.
   */
  private List<BedGraphGroupModel> mGroups = new ArrayList<BedGraphGroupModel>();

  /**
   * Add a BedGraph file to its own individual group.
   *
   * @param bedGraph the bed graph
   */
  public void add(UCSCTrack bedGraph) {
    add(DEFAULT_GROUP, bedGraph);
  }

  /**
   * Adds the.
   *
   * @param group    the group
   * @param bedGraph the bed graph
   */
  public void add(String group, UCSCTrack bedGraph) {
    createGroup(group);

    mGroups.get(mGroupMap.get(group)).add(bedGraph);
  }

  /**
   * Creates the group.
   *
   * @param group the group
   */
  public void createGroup(String group) {
    if (!mGroupMap.containsKey(group)) {
      BedGraphGroupModel model = new BedGraphGroupModel(group);

      model.addChangeListener(this);

      mGroupMap.put(group, mGroups.size());
      mGroups.add(model);
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.lib.event.ChangeListener#changed(org.abh.lib.event.ChangeEvent)
   */
  @Override
  public void changed(ChangeEvent e) {
    fireChanged();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<BedGraphGroupModel> iterator() {
    return mGroups.iterator();
  }

  /**
   * Gets the.
   *
   * @param group the group
   * @return the bed graph group model
   */
  public BedGraphGroupModel get(String group) {
    return mGroups.get(mGroupMap.get(group));
  }

  /**
   * Clear.
   */
  public void clear() {
    mGroups.clear();
    mGroupMap.clear();

    // fireChanged();
  }

  /**
   * Size.
   *
   * @return the int
   */
  public int size() {
    return mGroups.size();
  }

  /**
   * Removes the group.
   *
   * @param group the group
   */
  public void removeGroup(String group) {
    System.err.println("Remove group " + group);

    mGroups.remove((int) mGroupMap.get(group));
    mGroupMap.remove(group);

    fireChanged();
  }

  /**
   * Update group.
   *
   * @param oldGroup the old group
   * @param newGroup the new group
   */
  public void updateGroup(String oldGroup, String newGroup) {
    if (mGroupMap.containsKey(oldGroup)) {
      BedGraphGroupModel old = mGroups.get(mGroupMap.get(oldGroup));

      BedGraphGroupModel replacement = new BedGraphGroupModel(newGroup);

      for (UCSCTrack bedGraph : old) {
        replacement.add(bedGraph);
      }

      // update entries

      mGroups.add((int) mGroupMap.get(oldGroup), replacement);
      mGroups.remove((int) mGroupMap.get(oldGroup) + 1);

      mGroupMap.put(newGroup, mGroupMap.get(oldGroup));
      mGroupMap.remove(oldGroup);
    } else {
      createGroup(newGroup);
    }

    fireChanged();
  }
}