package org.jebtk.bioinformatics.genomic;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Mutation implements Comparable<Mutation> {
  private static final Pattern MUTATION_REGEX = Pattern.compile("([a-zA-Z]+)(\\d+)([a-zA-Z]+)");

  private final String mFrom;
  private final String mTo;
  private final MutationType mType;
  private final int mLocation;

  private int mIndex;

  public Mutation(int location, String from, String to) {
    mLocation = location;
    mFrom = from;
    mTo = to;
    mIndex = location - 1;

    if (from.length() > to.length()) {
      mType = MutationType.DELETION;
    } else if (from.length() < to.length()) {
      mType = MutationType.INSERTION;
    } else {
      mType = MutationType.CHANGE;
    }
  }

  public int getIndex() {
    return mIndex;
  }

  public int getLocation() {
    return mLocation;
  }

  public String getFrom() {
    return mFrom;
  }

  public String getTo() {
    return mTo;
  }

  public MutationType getType() {
    return mType;
  }

  @Override
  public String toString() {
    return mFrom + mLocation + mTo;
  }

  @Override
  public int compareTo(Mutation m) {
    if (mIndex > m.mIndex) {
      return 1;
    } else if (mIndex < m.mIndex) {
      return -1;
    } else {
      int c;

      c = mFrom.compareTo(m.mFrom);

      if (c != 0) {
        return c;
      } else {
        c = mTo.compareTo(m.mTo);

        if (c != 0) {
          return c;
        } else {
          return 0;
        }
      }
    }
  }

  @Override
  public boolean equals(Object o) {
    if (o instanceof Mutation) {
      return compareTo((Mutation) o) == 0;
    } else {
      return false;
    }
  }

  public static Mutation parse(String s) {
    Matcher matcher = MUTATION_REGEX.matcher(s);

    if (matcher.find()) {
      String from = matcher.group(1);
      int location = Integer.parseInt(matcher.group(2));
      String to = matcher.group(3);

      Mutation mutation = new Mutation(location, from, to);

      return mutation;
    } else {
      return null;
    }
  }

  public static List<Mutation> parse(Collection<String> lines) {
    List<Mutation> ret = new ArrayList<Mutation>(lines.size());

    for (String text : lines) {
      Mutation mutation = parse(text);

      if (mutation != null) {
        ret.add(mutation);
      }
    }

    return ret;
  }

  public static List<Mutation> parse(String[] lines) {
    List<Mutation> ret = new ArrayList<Mutation>(lines.length);

    for (String text : lines) {
      Mutation mutation = parse(text);

      if (mutation != null) {
        ret.add(mutation);
      }
    }

    return ret;
  }

  public static void toString(List<Mutation> mutations, PrintStream s) {
    for (Mutation mutation : mutations) {
      s.println(mutation.toString());
    }
  }

}
