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

import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.core.event.ChangeListeners;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.FormattedTxt;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

import com.fasterxml.jackson.annotation.JsonIgnore;

/**
 * Describes a region of a genome.
 *
 * @author Antony Holmes
 */
public class Region extends ChangeListeners implements Comparable<Region>, FormattedTxt {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /** The Constant REGION_REGEX. */
  private static final Pattern REGION_REGEX = Pattern.compile("\\d+(-\\d+)");

  public static final Region NULL_REGION = new Region();

  /**
   * The member start.
   */
  @JsonIgnore
  public final int mStart;

  /**
   * The member end.
   */
  @JsonIgnore
  public final int mEnd;

  /**
   * The member length.
   */
  @JsonIgnore
  public final int mLength;

  /**
   * Creates a special empty region
   */
  private Region() {
    mStart = 0;
    mEnd = 0;
    mLength = 0;
  }

  /**
   * Instantiates a new region.
   *
   * @param region the region
   */
  public Region(Region region) {
    this(region.mStart, region.mEnd);
  }

  /**
   * Instantiates a new region.
   *
   * @param start the start
   * @param end   the end
   */
  public Region(int start, int end) {
    // The start must be at least 1
    start = oneBased(start);
    end = oneBased(end);

    // Swap if the coordinates are the wrong way around
    if (start > end) {
      // int t = mStart;
      // mStart = mEnd;
      // mEnd = t;

      // hybrid
      start = start ^ end;
      // remove end to leave old start
      end = start ^ end;

      // start XOR start = end
      start = start ^ end;
    }

    mStart = start;
    // The end must be greater than the start
    mEnd = end; // Math.max(end, start + 1);
    mLength = mEnd - mStart + 1;
  }

  /**
   * Gets the start.
   *
   * @return the start
   */
  @JsonIgnore
  public int getStart() {
    return mStart;
  }

  /**
   * Gets the end.
   *
   * @return the end
   */
  @JsonIgnore
  public int getEnd() {
    return mEnd;
  }

  /**
   * Gets the length.
   *
   * @return the length
   */
  @JsonIgnore
  public int getLength() {
    return mLength;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  @JsonIgnore
  public String toString() {
    return toRange(mStart, mEnd);
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.lib.FormattedTxt#formattedTxt(java.lang.Appendable)
   */
  @Override
  public void formattedTxt(Appendable buffer) throws IOException {
    buffer.append(Integer.toString(mStart));
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Integer.toString(mEnd));
    buffer.append(TextUtils.NEW_LINE);
  }

  @Override
  public int compareTo(Region r) {
    if (mStart > r.mStart) {
      return 1;
    } else if (mStart < r.mStart) {
      return -1;
    } else {
      // If the starts are the same, rank by end coordinate
      if (mEnd > r.mEnd) {
        return 1;
      } else if (mEnd < r.mEnd) {
        return -1;
      } else {
        // both start and end equal each other
        return 0;
      }
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof Region) {
      return compareTo((Region) o) == 0;
    } else {
      return false;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    return toString().hashCode();
  }

  /**
   * Returns true if the text matches a region of the form m-n.
   *
   * @param region the region
   * @return true, if is region
   */
  protected static boolean isRegion(String region) {
    return REGION_REGEX.matcher(region).matches();
  }

  public static String region(String text) {
    if (text == null) {
      return TextUtils.EMPTY_STRING;
    }

    Matcher matcher = REGION_REGEX.matcher(text);

    if (matcher.find()) {
      return matcher.group(1);
    } else {
      return TextUtils.EMPTY_STRING;
    }
  }

  public static Region parseRegion(String location) {
    if (Io.isEmptyLine(location)) {
      return null;
    }

    if (location.contains(TextUtils.NA)) {
      return null;
    }

    if (location.length() == 0) {
      return null;
    }

    if (isRegion(location)) {
      location = region(location);

      int start;
      int end;

      if (location.indexOf("-") != -1) {
        List<String> tokens = Splitter.on('-').text(location); // )
                                                               // .(tokens.get(1),
                                                               // '-');

        start = TextUtils.parseInt(tokens.get(0));
        end = TextUtils.parseInt(tokens.get(1));
      } else {
        // single position

        start = TextUtils.parseInt(location);
        end = start;
      }

      return createRegion(start, end);
    } else {
      return null;
    }
  }

  public static Region createRegion(int start, int end) {
    return new Region(start, end);
  }

  /**
   * Shift.
   *
   * @param region the region
   * @param shift  the shift
   * @return the region
   */
  public static Region shift(Region region, int shift) {
    // bound the positions so they dont exceed the chromosome bounds
    int start = oneBased(region.mStart + shift);
    int end = oneBased(region.mEnd + shift);

    return new Region(start, end);
  }

  public static String toRange(int start, int end) {
    return start + "-" + end;
  }

  /**
   * Simply ensure coordinate is a minimum of one
   * 
   * @param v
   * @return
   */
  public static int oneBased(int v) {
    return Math.max(1, v);
  }

  /**
   * Ensure that end is
   * 
   * @param start
   * @param end
   * @return
   */
  public static int strictEnd(int start, int end) {
    return Math.max(start + 1, end);
  }
}
