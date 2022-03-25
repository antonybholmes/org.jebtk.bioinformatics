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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.search.Feature;
import org.jebtk.core.TableData;

/**
 * The class Search.
 */
public class Search {

  /**
   * Instantiates a new search.
   */
  private Search() {
    // do nothing
  }

  /**
   * Performs a binary search on an ordered Feature list allowing for precise or
   * range matching.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static final int binaryRangeSearch(int value, List<Feature> features) {
    int min = 0;
    int max = features.size() - 1;

    int mid = -1;

    // int midMax = features.size() - 1;

    while (min <= max) {
      mid = (min + max) / 2;

      // System.out.println("value:" + value + "Feature: " +
      // features.get(mid).getStart() + " " + features.get(mid).getEnd() + " " +
      // features.get(mid).getName());

      if (features.get(mid).getStart() <= value && features.get(mid).getEnd() >= value) {
        return mid;
      } else if (value > features.get(mid).getEnd()) {
        min = mid + 1;
      } else {
        max = mid - 1;
      }
    }

    return -1;
  }

  /**
   * Binary range search inner.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static final int binaryRangeSearchInner(int value, List<Feature> features) {
    int min = 0;
    int max = features.size() - 1;

    int mid = -1;

    int midMax = features.size() - 1;

    while (min <= max) {
      mid = (min + max) / 2;

      // System.out.println("value:" + value + "Feature: " +
      // features.get(mid).getStart() + " " + features.get(mid).getEnd() + " " +
      // features.get(mid).getName());

      if (features.get(mid).getStart() <= value && features.get(mid).getEnd() >= value) {
        return mid;
      } else if (mid < midMax && features.get(mid).getEnd() < value && features.get(mid + 1).getStart() > value) {
        // nearest position lies between two features
        // so pick the inner one (inclusive) {
        return mid + 1;
      } else if (value > features.get(mid).getEnd()) {
        min = mid + 1;
      } else {
        max = mid - 1;
      }
    }

    return -1;
  }

  /**
   * Binary range search outer.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static final int binaryRangeSearchOuter(int value, List<Feature> features) {
    int min = 0;
    int max = features.size() - 1;

    int mid = -1;

    int midMax = features.size() - 1;

    while (min <= max) {
      mid = (min + max) / 2;

      // System.out.println("value:" + value + "Feature: " +
      // features.get(mid).getStart() + " " + features.get(mid).getEnd() + " " +
      // features.get(mid).getName());

      if (features.get(mid).getStart() <= value && features.get(mid).getEnd() >= value) {
        return mid;
      } else if (mid < midMax && features.get(mid).getEnd() < value && features.get(mid + 1).getStart() > value) {
        return mid;
      } else if (value > features.get(mid).getEnd()) {
        min = mid + 1;
      } else {
        max = mid - 1;
      }
    }

    return -1;
  }

  /**
   * Skips over blocks of data until the data is passed and then moves back a
   * block and a more thorough search is performed. Required for genes etc where
   * they are ordered but overlapping so a binary search will not work
   *
   * @param value    the value
   * @param features the features
   * @param skip     the skip
   * @return the int
   */
  public static final int jumpSearchInner(int value, List<Feature> features, int skip) {
    int index = -1;

    boolean stop = false;
    boolean last = false;

    int max = -1;

    int i1 = 0;

    while (!stop) {
      // System.out.println(features.get(i).getName() + " " +
      // features.get(i).getStart() + " " + value);

      if (features.get(i1).getStart() > value) {
        // reached a point beyond what we want so trace back and search more
        // thoroughly

        // System.out.println("test" + features.get(i).getName() + " " +
        // features.get(i).getStart() );

        // go back
        int start = Math.max(0, i1 - skip);

        for (int j = start; j <= i1; ++j) {
          // System.out.println("real" + features.get(j).getName() + " " +
          // features.get(j).getStart() + " " + max);
          // find the greatest start of a segment spanning this point
          if (features.get(j).getStart() <= value && features.get(j).getEnd() >= value
              && features.get(j).getStart() >= max) {
            index = j;

            max = features.get(j).getStart();
          }
        }

        if (index == -1) {
          // there was nothing spanning so pick the closest inside
          // the range
          for (int j = start; j <= i1; ++j) {
            if (features.get(j).getStart() > value) {
              index = j;
              break;
            }
          }
        }

        break;
      }

      if (last) {
        stop = true;
        break;
      }

      i1 += skip;

      if (i1 > features.size()) {
        i1 = features.size() - 1;
        last = true;
      }
    }

    // System.out.println("returning " + index);

    return index;
  }

  /**
   * Jump search outer.
   *
   * @param value    the value
   * @param features the features
   * @param skip     the skip
   * @return the int
   */
  public static final int jumpSearchOuter(int value, List<Feature> features, int skip) {
    int index = -1;

    int max = -1;

    boolean stop = false;
    boolean last = false;

    int i = 0;

    while (!stop) {
      if (features.get(i).getStart() > value) {
        // reached a point beyond what we want so trace back and search more
        // thoroughly

        // System.out.println(features.get(i).getName() + " " +
        // features.get(i).getStart() + " " + value);

        int start = Math.max(0, i - skip);

        for (int j = start; j <= i; ++j) {
          if (features.get(j).getStart() <= value && features.get(j).getEnd() >= max) {
            index = j;

            max = features.get(j).getEnd();
          }
        }

        break;
      }

      if (last) {
        stop = true;
        break;
      }

      i += skip;

      if (i > features.size()) {
        i = features.size() - 1;
        last = true;
      }
    }

    return index;
  }

  /**
   * Performs a basic linear search of the features and return the first position
   * to make the cut.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static int find(int value, List<Feature> features) {
    for (int i = 0; i < features.size(); ++i) {
      if (features.get(i).getStart() <= value && features.get(i).getEnd() >= value) {
        return i;
      }
    }

    return -1;
  }

  /**
   * Find inner.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static int findInner(int value, List<Feature> features) {
    int min = Integer.MAX_VALUE;
    int mini = -1;

    for (int i = 0; i < features.size(); ++i) {
      if (features.get(i).getStart() <= value && features.get(i).getEnd() >= value) {
        return i;
      }

      if (features.get(i).getStart() >= value) {
        if (features.get(i).getStart() - value < min) {
          min = features.get(i).getStart() - value;
          mini = i;
        }
      }
    }

    return mini;
  }

  /**
   * Find outer.
   *
   * @param value    the value
   * @param features the features
   * @return the int
   */
  public static int findOuter(int value, List<Feature> features) {
    int min = Integer.MAX_VALUE;
    int mini = -1;

    for (int i = 0; i < features.size(); ++i) {
      if (features.get(i).getStart() <= value && features.get(i).getEnd() >= value) {
        return i;
      }

      if (features.get(i).getStart() <= value) {
        if (value - features.get(i).getStart() < min) {
          min = value - features.get(i).getStart();
          mini = i;
        }
      }
    }

    return mini;
  }

  /**
   * Creates a Feature list with single position features (i.e. the start and end
   * points are the same on a given Feature.
   *
   * @param table         the table
   * @param startColumn   the start column
   * @param groupColumn   the group column
   * @param featureColumn the feature column
   * @return the map
   */
  public static final Map<Integer, List<Feature>> createFeatureList(TableData<String> table, int startColumn,
      int groupColumn, int featureColumn) {
    return createFeatureList(table, startColumn, startColumn, groupColumn, featureColumn);
  }

  /**
   * Creates a Feature map from a table.
   *
   * @param table         a table to create Feature data from
   * @param startColumn   the start positions
   * @param endColumn     the end positions
   * @param groupColumn   features can be grouped, such as by chromosome
   * @param featureColumn identifier such as a gene name
   * @return the map
   */
  public static final Map<Integer, List<Feature>> createFeatureList(TableData<String> table, int startColumn,
      int endColumn, int groupColumn, int featureColumn) {
    int group;

    Map<Integer, List<Feature>> features = new HashMap<Integer, List<Feature>>();

    for (List<String> row : table) {
      try {
        group = Integer.parseInt(row.get(groupColumn));

        if (!features.containsKey(group)) {
          features.put(group, new ArrayList<Feature>());
        }

        Feature feature = new Feature(featureColumn != -1 ? row.get(featureColumn) : null, Chromosome.NO_CHR,
            Integer.parseInt(row.get(startColumn)), Integer.parseInt(row.get(endColumn)));

        // if (group == 13) {
        // System.out.println("Feature: " + Feature.start + " " + Feature.end +
        // " " +
        // Feature.type);
        // }

        features.get(group).add(feature);
      } catch (Exception e) {
        // e.printStackTrace();
      }
    }

    return features;
  }

  /*
   * public static final Set<String> features(Map<Integer, List<Feature>>
   * allfeatures, int start, int end, int group) { int indexStart; int indexEnd;
   * 
   * List<Feature> features = allfeatures.get(group);
   * 
   * Set<String> ret = new HashSet<String>();
   * 
   * if (features == null) { return null; }
   * 
   * indexStart = Search.binaryRangeSearch(start, features, SearchType.INNER);
   * 
   * if (indexStart == -1) { indexStart = 0; }
   * 
   * // TODO what is the correct place to return an index from indexEnd =
   * Search.binaryRangeSearch(end, features, SearchType.OUTER);
   * 
   * if (indexEnd == -1) { indexEnd = features.size() - 1; }
   * 
   * if (indexEnd < indexStart) { indexEnd = indexStart; }
   * 
   * for (int j = indexStart; j <= indexEnd; ++j) {
   * //LogServer.getInstance().getLogger(getApplicationInformation().getName()).
   * write("s: " + snpIndexStart + " e:" + snpIndexEnd + " " + ratio + " " +
   * chromosomefeatures.get(j).type); ret.add(features.get(j).getName()); }
   * 
   * return ret; }
   */
}
