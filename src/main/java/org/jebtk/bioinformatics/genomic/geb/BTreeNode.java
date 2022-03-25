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
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

/**
 * Represents a node in a radix tree. Names are stored as sequences of
 * characters which can be associated with one or more objects.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class BTreeNode<T> implements Comparable<BTreeNode<T>>, Serializable, Iterable<T> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The member objects.
   */
  // private List<T> mObjects = new UniqueArrayList<T>();

  private Set<T> mObjects = new TreeSet<T>();

  /**
   * The member children.
   */
  private BTreeNode<T> mC1;
  private BTreeNode<T> mC2;

  /**
   * The member c.
   */
  private int mBin;

  /**
   * Instantiates a new radix object node.
   *
   * @param c      the c
   * @param prefix the prefix
   */
  public BTreeNode(int bin) {
    mBin = bin;
  }

  /**
   * Gets the char.
   *
   * @return the char
   */
  public int getIndex() {
    return mBin;
  }

  public void setC1(BTreeNode<T> c) {
    mC1 = c;
  }

  public void setC2(BTreeNode<T> c) {
    mC2 = c;
  }

  public BTreeNode<T> getC1() {
    return mC1;
  }

  public BTreeNode<T> getC2() {
    return mC2;
  }

  /**
   * Gets the objects.
   *
   * @return the objects
   */
  public Iterable<T> getObjects() {
    return mObjects;
  }

  /**
   * Returns the node associated with a given prefix.
   *
   * @param prefix the prefix
   * @return the child
   */
  public BTreeNode<T> getChild(int bin) {
    return getChild(this, bin);
  }

  private static <TT> BTreeNode<TT> getChild(BTreeNode<TT> root, int bin) {
    while (root != null) {
      int b = root.getIndex();

      if (b == bin) {
        return root;
      } else if (bin < b) {
        root = root.getC1();
      } else {
        root = root.getC2();
      }
    }

    return null;
  }

  public void add(T v) {
    add(mBin, v);
  }

  /**
   * Parse a string into prefixs and build a sub tree under the current node to
   * represent that string.
   *
   * @param word the word
   * @param v    the object
   */
  public void add(int bin, T v) {
    if (bin == mBin) {
      mObjects.add(v);
    } else {
      BTreeNode<T> ret = new BTreeNode<T>(bin);

      if (bin < mBin) {
        setC1(ret);
      } else {
        setC2(ret);
      }
    }
  }

  public int getObjectCount() {
    return mObjects.size();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(BTreeNode<T> n) {
    if (mBin == n.mBin) {
      return 0;
    } else if (mBin < n.mBin) {
      return -1;
    } else {
      return 1;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof BTreeNode)) {
      return false;
    }

    return compareTo((BTreeNode<T>) o) == 0;
  }

  /**
   * Clear.
   */
  public void clear() {
    mObjects.clear();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<T> iterator() {
    return mObjects.iterator();
  }
}
