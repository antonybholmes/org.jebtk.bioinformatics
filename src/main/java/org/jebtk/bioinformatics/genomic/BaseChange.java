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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.core.text.TextUtils;

/**
 * Describes a base change at a position.
 *
 * @author Antony Holmes
 *
 */
public class BaseChange implements Comparable<BaseChange> {

  /**
   * The constant BASE_CHANGE_REGEX_TEXT.
   */
  public static final String BASE_CHANGE_REGEX_TEXT = "^([^\\d]+?)(\\d+)([^\\d]+)";

  /**
   * The constant BASE_CHANGE_REGEX.
   */
  public static final Pattern BASE_CHANGE_REGEX = Pattern.compile(BASE_CHANGE_REGEX_TEXT);

  /** The Constant MUTATION_REGEX_TEXT. */
  public static final String MUTATION_REGEX_TEXT = ".*([A-Z])(\\d+)([A-Z]).*";

  /** The Constant MUTATION_REGEX. */
  public static final Pattern MUTATION_REGEX = Pattern.compile(MUTATION_REGEX_TEXT);

  /** The Constant DELETION_REGEX_TEXT. */
  public static final String DELETION_REGEX_TEXT = ".*([A-Z])(\\d+)_([A-Z])(\\d+)del.*";

  /** The Constant DELETION_REGEX. */
  public static final Pattern DELETION_REGEX = Pattern.compile(DELETION_REGEX_TEXT);

  /**
   * The constant RANGE_CHANGE_REGEX_TEXT.
   */
  public static final String RANGE_CHANGE_REGEX_TEXT = "^([^\\d]+)\\((\\d+)-(\\d+)\\)([^\\d]+)";

  /**
   * The constant RANGE_CHANGE_REGEX.
   */
  public static final Pattern RANGE_CHANGE_REGEX = Pattern.compile(RANGE_CHANGE_REGEX_TEXT);

  /**
   * The from.
   */
  private String mFrom;

  /**
   * The location.
   */
  private int mLocation;

  /**
   * The to.
   */
  private String mTo;

  /**
   * The text.
   */
  private String mText;

  /**
   * Instantiates a new base change.
   *
   * @param text the text
   */
  public BaseChange(String text) {
    text = TextUtils.chomp(text);

    if (text.matches(RANGE_CHANGE_REGEX_TEXT)) {
      Matcher matcher = RANGE_CHANGE_REGEX.matcher(text);

      matcher.find();

      int l = (Integer.parseInt(matcher.group(2)) + Integer.parseInt(matcher.group(3))) / 2;

      setup(matcher.group(1), matcher.group(4), l);

    } else if (text.matches(BASE_CHANGE_REGEX_TEXT)) {
      Matcher matcher = BASE_CHANGE_REGEX.matcher(text);

      matcher.find();

      setup(matcher.group(1), matcher.group(3), Integer.parseInt(matcher.group(2)));
    } else if (text.matches(MUTATION_REGEX_TEXT)) {
      Matcher matcher = MUTATION_REGEX.matcher(text);

      matcher.find();

      setup(matcher.group(1), matcher.group(3), Integer.parseInt(matcher.group(2)));
    } else if (text.matches(DELETION_REGEX_TEXT)) {
      Matcher matcher = DELETION_REGEX.matcher(text);

      matcher.find();

      setup(matcher.group(1), "-", Integer.parseInt(matcher.group(2)));
    } else {
      setup(text, text, -1);
    }
  }

  /**
   * Create a base change.
   *
   * @param from     the from
   * @param to       the to
   * @param location the location
   */
  public BaseChange(String from, String to, int location) {
    mFrom = from;
    mTo = to;
    mLocation = location;

    setup(from, to, location);
  }

  /**
   * Setup.
   *
   * @param from     the from
   * @param to       the to
   * @param location the location
   */
  private void setup(String from, String to, int location) {
    this.mFrom = from;
    this.mTo = to;
    this.mLocation = location;

    this.mText = from + Integer.toString(location) + to;
  }

  /**
   * Gets the from.
   *
   * @return the from
   */
  public final String getFrom() {
    return mFrom;
  }

  /**
   * Gets the to.
   *
   * @return the to
   */
  public final String getTo() {
    return mTo;
  }

  /**
   * Gets the location.
   *
   * @return the location
   */
  public final int getLocation() {
    return mLocation;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  public String toString() {
    return mText;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  public int compareTo(BaseChange b) {
    return mText.compareTo(b.toString());
  }
}
