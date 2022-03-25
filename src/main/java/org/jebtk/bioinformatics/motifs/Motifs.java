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

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.jebtk.core.NameGetter;
import org.jebtk.core.collections.CollectionUtils;

/**
 * Represents a motif sequence. Each position consists of the counts for each
 * base (a,c,g,t).
 * 
 * @author Antony Holmes
 *
 */
public class Motifs implements NameGetter, Iterable<Motif> {

  /** The m name. */
  private String mName;

  /** The m motifs. */
  private List<Motif> mMotifs;

  /**
   * Instantiates a new motifs.
   *
   * @param name   the name
   * @param motifs the motifs
   */
  public Motifs(String name, Collection<Motif> motifs) {
    mName = name;

    mMotifs = CollectionUtils.sort(motifs);
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.common.NameProperty#getName()
   */
  @Override
  public String getName() {
    return mName;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<Motif> iterator() {
    return mMotifs.iterator();
  }
}
