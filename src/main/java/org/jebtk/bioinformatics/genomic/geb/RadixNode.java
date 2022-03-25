/**
 * Copyright 2016 Antony Holmes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.jebtk.bioinformatics.genomic.geb;

import java.io.Serializable;
import java.util.Collection;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;

/**
 * Represents a node in a radix tree. Names are stored as sequences of
 * characters which can be associated with one or more objects.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class RadixNode<T> implements Comparable<RadixNode<T>>, Serializable {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The member objects.
   */
  // private List<T> mObjects = new UniqueArrayList<T>();

  private Set<T> mExactObjects = new TreeSet<T>();
  private Set<T> mObjects = new TreeSet<T>();

  /**
   * The member children.
   */
  private IterMap<Character, RadixNode<T>> mChildren = new IterTreeMap<Character, RadixNode<T>>();

  /**
   * The member c.
   */
  private char mC;

  /**
   * The member prefix.
   */
  // private String mPrefix;

  /**
   * Create the root node
   */
  public RadixNode() {
    this((char) 0);
  }

  /**
   * Instantiates a new radix object node.
   *
   * @param c      the c
   * @param prefix the prefix
   */
  public RadixNode(char c) {
    mC = standardize(c);
    // mPrefix = newPrefix(prefix, c);
  }

  /**
   * Gets the char.
   *
   * @return the char
   */
  public char getChar() {
    return mC;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  public String toString() {
    return Character.toString(mC);
  }

  /**
   * Gets the objects.
   *
   * @return the objects
   */
  public Collection<T> getExactObjects() {
    return mExactObjects;
  }

  public Collection<T> getObjects() {
    return mObjects;
  }

  /**
   * Returns the node associated with a given prefix.
   *
   * @param prefix the prefix
   * @return the child
   */
  public RadixNode<T> getChild(String prefix) {
    return getChild(this, prefix);
  }

  private static <TT> RadixNode<TT> getChild(RadixNode<TT> root, String prefix) {
    RadixNode<TT> ret = root;

    char[] chars = standardize(prefix).toCharArray();

    for (char c : chars) {
      // Create child
      if (!ret.mChildren.containsKey(c)) {
        ret.mChildren.put(c, new RadixNode<TT>(c));
      }

      ret = ret.mChildren.get(c);
    }

    return ret;
  }

  /**
   * Returns a child node. Will return null if the child does not exist.
   *
   * @param c the c
   * @return the child
   */
  public RadixNode<T> getChild(char c) {
    return mChildren.get(standardize(c));
  }

  /**
   * Parse a string into prefixs and build a sub tree under the current node to
   * represent that string.
   *
   * @param word the word
   * @param v    the object
   */
  public void add(String word, T v) {
    // Add object to this node
    mObjects.add(v);

    RadixNode<T> ret = this;

    char[] chars = standardize(word).toCharArray();

    for (int i = 0; i < chars.length; ++i) {
      char c = chars[i];

      // Create child
      if (!ret.mChildren.containsKey(c)) {
        ret.mChildren.put(c, new RadixNode<T>(c));
      }

      ret = ret.mChildren.get(c);

      if (i == chars.length - 1) {
        ret.mExactObjects.add(v);
      } else {
        ret.mObjects.add(v);
      }
    }
  }

  /**
   * Gets the child count.
   *
   * @return the child count
   */
  public int getChildCount() {
    return mChildren.size();
  }

  public Set<Entry<Character, RadixNode<T>>> getChildren() {
    return mChildren.entrySet();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(RadixNode<T> n) {
    if (mC > n.mC) {
      return 1;
    } else if (mC < n.mC) {
      return -1;
    } else {
      return 0;
    }
  }

  /**
   * Clear.
   */
  public void clear() {
    mChildren.clear();

    mObjects.clear();
  }

  /**
   * Ensure characters are consistent for searching purposes i.e case insensitive.
   *
   * @param name the name
   * @return the char
   */
  private static char standardize(char name) {
    return Character.toLowerCase(name);
  }

  /**
   * Standardize.
   *
   * @param name the name
   * @return the string
   */
  private static String standardize(String name) {
    return name.toLowerCase();
  }

  public Set<Character> getChildNames() {
    return mChildren.keySet();
  }
}
