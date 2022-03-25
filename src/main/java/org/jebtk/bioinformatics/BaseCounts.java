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

import java.util.Arrays;

import org.jebtk.core.Mathematics;

/**
 * Stores the counts of each base (a,c,g,t).
 * 
 * @author Antony Holmes
 *
 */
public class BaseCounts {

  /**
   * The member counts.
   */
  private double[] mCounts = new double[5];
  private double mMax;

  /**
   * Instantiates a new base counts.
   *
   * @param a the a
   * @param c the c
   * @param g the g
   * @param t the t
   */
  public BaseCounts(double a, double c, double g, double t) {
    this(a, c, g, t, 0);
  }

  public BaseCounts(double a, double c, double g, double t, double n) {
    this(a, c, g, t, n, true);
  }

  public BaseCounts(double a, double c, double g, double t, boolean norm) {
    this(a, c, g, t, 0, norm);
  }

  public BaseCounts(double a, double c, double g, double t, double n, boolean norm) {

    mCounts[0] = a;
    mCounts[1] = c;
    mCounts[2] = g;
    mCounts[3] = t;
    mCounts[4] = n;

    if (norm) {
      double total = a + c + g + t + n;

      mCounts[0] /= total;
      mCounts[1] /= total;
      mCounts[2] /= total;
      mCounts[3] /= total;
      mCounts[4] /= total;
    }

    mMax = Mathematics.max(mCounts[0], mCounts[1], mCounts[2], mCounts[3], mCounts[4]);

  }

  /**
   * Gets the count.
   *
   * @param base the base
   * @return the count
   */
  public double getCount(char base) {
    switch (base) {
    case 'A':
    case 'a':
      return getCount(0);
    case 'C':
    case 'c':
      return getCount(1);
    case 'G':
    case 'g':
      return getCount(2);
    case 'T':
    case 't':
    case 'U':
    case 'u':
      return getCount(3);
    default:
      // N etc
      return getCount(4);// getCount(0); //0.25;
    }
  }

  /**
   * Gets the count.
   *
   * @param base the base
   * @return the count
   */
  public double getCount(int base) {
    return mCounts[base];
  }

  /**
   * Gets the a.
   *
   * @return the a
   */
  public double getA() {
    return getCount(0);
  }

  /**
   * Gets the c.
   *
   * @return the c
   */
  public double getC() {
    return getCount(1);
  }

  /**
   * Gets the g.
   *
   * @return the g
   */
  public double getG() {
    return getCount(2);
  }

  /**
   * Gets the t.
   *
   * @return the t
   */
  public double getT() {
    return getCount(3);
  }

  public double getU() {
    return getCount(3);
  }

  public double getN() {
    return getCount(4);
  }

  public double getMaxScore() {
    return mMax;
  }

  @Override
  public String toString() {
    return Arrays.toString(mCounts);
  }

  /**
   * Returns the sum of all bases.
   * 
   * @return the sum of all bases.
   */
  public double getSum() {
    return Mathematics.sum(mCounts);
  }
}
