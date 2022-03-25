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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * The class AminoAcid.
 */
public class AminoAcid implements Comparable<AminoAcid>, Iterable<Codon> {

  /**
   * The abbreviation.
   */
  private String abbreviation;

  /**
   * The name.
   */
  private String name;

  /**
   * The codons.
   */
  private List<Codon> codons = new ArrayList<Codon>();

  /**
   * The display.
   */
  private String display;

  /**
   * The letter.
   */
  private char letter;

  /**
   * Instantiates a new amino acid.
   *
   * @param name         the name
   * @param abbreviation the abbreviation
   * @param letter       the letter
   */
  public AminoAcid(String name, String abbreviation, char letter) {
    this.name = name;
    this.abbreviation = abbreviation;
    this.letter = letter;

    this.display = name + "/" + abbreviation + "/" + letter;
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public String getName() {
    return name;
  }

  /**
   * Gets the abbreviation.
   *
   * @return the abbreviation
   */
  public String getAbbreviation() {
    return abbreviation;
  }

  /**
   * Gets the letter.
   *
   * @return the letter
   */
  public char getLetter() {
    return letter;
  }

  /**
   * Adds the codon.
   *
   * @param codon the codon
   */
  public void addCodon(Codon codon) {
    codons.add(codon);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return display;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<Codon> iterator() {
    return codons.iterator();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(AminoAcid a) {
    return name.compareTo(a.name);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof AminoAcid)) {
      return false;
    }

    return compareTo((AminoAcid) o) == 0;
  }
}
