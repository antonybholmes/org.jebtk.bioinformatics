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
package org.jebtk.bioinformatics.genomic;

/**
 * Represents a gene along with a description of the gene.
 * 
 * @author Antony Holmes
 *
 */
public class ProbeGene implements Comparable<ProbeGene> {

  /**
   * The member name.
   */
  private String mName;

  /**
   * The member description.
   */
  private String mGeneSymbol;

  private String mDescription;

  /**
   * Instantiates a new probe gene.
   *
   * @param name        the name
   * @param description the description
   */
  public ProbeGene(String name, String symbol) {
    this(name, symbol, symbol);
  }

  /**
   * Instantiates a new probe gene.
   *
   * @param name        the name
   * @param description the description
   * @param chr         the chr
   */
  public ProbeGene(String name, String symbol, String description) {
    mName = name;
    mGeneSymbol = symbol;
    mDescription = description;
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public String getName() {
    return mName;
  }

  /**
   * Gets the gene symbol.
   *
   * @return the gene symbol.
   */
  public String getGeneSymbol() {
    return mGeneSymbol;
  }

  /**
   * Gets the description.
   *
   * @return the description
   */
  public String getDescription() {
    return mDescription;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return mName;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(ProbeGene gene) {
    return mName.compareTo(gene.mName);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof ProbeGene)) {
      return false;
    }

    return compareTo((ProbeGene) o) == 0;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    return mName.hashCode();
  }
}
