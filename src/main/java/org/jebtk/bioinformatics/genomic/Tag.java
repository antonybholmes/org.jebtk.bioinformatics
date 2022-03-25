package org.jebtk.bioinformatics.genomic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jebtk.core.text.TextUtils;

import com.fasterxml.jackson.annotation.JsonGetter;
import com.fasterxml.jackson.annotation.JsonIgnore;

public abstract class Tag implements Comparable<Tag> {
  @JsonIgnore
  public int toInt() {
    return Integer.MIN_VALUE;
  }

  @JsonIgnore
  public double toDouble() {
    return Double.MIN_VALUE;
  }

  @JsonIgnore
  public abstract TagType getType();

  @Override
  @JsonGetter("value")
  public String toString() {
    return super.toString();
  }

  @Override
  public int compareTo(Tag p) {
    int ret = getType().compareTo(p.getType());

    if (ret != 0) {
      return ret;
    }

    return compareValue(p);
  }

  /**
   * Order Props by their value if the names are not sufficient.
   * 
   * @param p
   * @return
   */
  protected abstract int compareValue(Tag p);

  /**
   * Standardize property names.
   * 
   * @param name
   * @return
   */
  public static String format(String name) {
    return name.toLowerCase();
  }

  public static List<Tag> toTags(Collection<String> tags) {
    ArrayList<Tag> ret = new ArrayList<Tag>();

    for (String tag : tags) {
      ret.add(toTag(tag));
    }

    return ret;
  }

  private static Tag toTag(String tag) {
    if (TextUtils.isInt(tag)) {
      return new IntTag(Integer.parseInt(tag));
    } else if (TextUtils.isDouble(tag)) {
      return new IntTag(Integer.parseInt(tag));
    } else {
      return new TextTag(tag);
    }
  }
}
