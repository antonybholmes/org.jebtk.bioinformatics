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

import org.jebtk.core.event.ChangeEvent;
import org.jebtk.core.event.ChangeListener;
import org.jebtk.core.model.ListModel;

/**
 * A model for shared bed graphs.
 * 
 * @author Antony Holmes
 *
 */
public class BedGraphGroupModel extends ListModel<UCSCTrack> implements ChangeListener {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The member name.
   */
  private String mName;

  /**
   * Instantiates a new bed graph group model.
   *
   * @param name the name
   */
  public BedGraphGroupModel(String name) {
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
   * @see org.abh.lib.model.ListModel#add(java.lang.Object)
   */
  @Override
  public boolean add(UCSCTrack bedGraph) {
    boolean ret = super.add(bedGraph);

    bedGraph.addChangeListener(this);

    return ret;
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
}