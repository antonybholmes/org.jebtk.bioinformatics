/**
 * Copyright (c) 2016, Antony Holmes
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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Describes a genomic locations strand.
 */
public enum Strand {
  /// ** Unspecified */
  // NONE(0),

  /** The sense. */
  SENSE(0),

  /** The antisense. */
  ANTISENSE(1);

  /** The strand's value. */
  private int mValue;

  /**
   * Instantiates a new strand.
   *
   * @param value the value
   */
  private Strand(int value) {
    mValue = value;
  }

  /**
   * Gets the value.
   *
   * @return the value
   */
  public int getValue() {
    return mValue;
  }

  @Override
  public String toString() {
    return toString(this);
  }

  /**
   * Parses the.
   *
   * @param strand the strand
   * @return the strand
   */
  public static Strand parse(char strand) {
    switch (strand) {
    case '-':
      return ANTISENSE;
    default:
      return SENSE;
    }
  }

  /**
   * Parses the.
   *
   * @param strand the strand
   * @return the strand
   */
  public static Strand parse(String strand) {
    return parse(strand.charAt(0));
  }

  /**
   * Convert a strand to its character equivalent: SENSE -> '+' ANTI_SENSE -> '-'
   * NONE -> '.'
   *
   * @param strand the strand
   * @return the string
   */
  public static String toString(Strand strand) {
    return Character.toString(toChar(strand));
  }

  /**
   * To char.
   *
   * @param strand the strand
   * @return the char
   */
  public static char toChar(Strand strand) {
    switch (strand) {
    case ANTISENSE:
      return '-';
    default:
      return '+';
    }
  }

  /**
   * Convert a collection of strands to their character equivalents.
   *
   * @param strands the strands
   * @return the list
   */
  public static List<Character> toChar(final Collection<Strand> strands) {
    List<Character> ret = new ArrayList<Character>(strands.size());

    for (Strand s : strands) {
      ret.add(toChar(s));
    }

    return ret;
  }

  /**
   * Parse a list of characters as strands.
   *
   * @param chars the chars
   * @return the list
   */
  public static List<Strand> parse(final Collection<Character> chars) {
    List<Strand> ret = new ArrayList<Strand>(chars.size());

    for (char c : chars) {
      ret.add(parse(c));
    }

    return ret;
  }

  public static Strand[] parse(final char[] chars) {
    Strand[] ret = new Strand[chars.length];

    for (int i = 0; i < chars.length; ++i) {
      ret[i] = parse(chars[i]);
    }

    return ret;
  }

  /**
   * Checks if is sense.
   *
   * @param strand the strand
   * @return true, if is sense
   */
  public static boolean isSense(Strand strand) {
    return strand == SENSE;
  }

  /**
   * Return the opposite strand.
   * 
   * @param s
   * @return
   */
  public static Strand oppositeStrand(Strand s) {
    switch (s) {
    case ANTISENSE:
      return SENSE;
    default:
      return ANTISENSE;
    }
  }
}
