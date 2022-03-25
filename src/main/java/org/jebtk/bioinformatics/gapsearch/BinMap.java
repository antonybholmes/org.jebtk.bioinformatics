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
package org.jebtk.bioinformatics.gapsearch;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.jebtk.core.Range;
import org.jebtk.core.collections.IterMap;

/**
 * Use fixed size blocks to find features.
 *
 * @author Antony Holmes
 * @param <T> the generic type
 */
public class BinMap<T> implements IterMap<Integer, T> {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant DEFAULT_BIN_SIZE.
   */
  public static final int DEFAULT_BIN_SIZE = 10000;

  /**
   * A chromosome does not have more than 270 million bases
   */
  public static final int MAX_BP = 270000000;

  /**
   * The member bin size.
   */
  protected final int mBinSize;

  private Object[] mData;

  private class BinMapIter implements Iterator<Entry<Integer, T>> {

    private int mC = 0;

    @Override
    public boolean hasNext() {
      return mC < mData.length;
    }

    @Override
    public Entry<Integer, T> next() {
      Entry<Integer, T> e = new org.jebtk.core.collections.Entry<Integer, T>(mC, (T) mData[mC]);
      ++mC;

      return e;
    }

    @Override
    public void remove() {
      // TODO Auto-generated method stub

    }

  }

  /**
   * Instantiates a new fixed gap search.
   */
  public BinMap() {
    this(DEFAULT_BIN_SIZE);
  }

  /**
   * Instantiates a new fixed gap search.
   *
   * @param binSize the bin size
   */
  public BinMap(int binSize) {
    mBinSize = Math.max(1, binSize);
    mData = new Object[MAX_BP / mBinSize];
  }

  @Override
  public void clear() {
    Arrays.fill(mData, null);
  }

  @Override
  public boolean containsKey(Object o) {
    if (o instanceof Integer) {
      return containsKey((Integer) o);
    } else {
      return false;
    }
  }

  public boolean containsKey(Integer start) {
    return containsKey((int) start);
  }

  public boolean containsKey(int start) {
    return mData[start] != null;
  }

  @Override
  public boolean containsValue(Object arg0) {
    return false;
  }

  @Override
  public Set<java.util.Map.Entry<Integer, T>> entrySet() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public T get(Object o) {
    if (o instanceof Integer) {
      return get((Integer) o);
    } else {
      return null;
    }
  }

  public T get(Integer start) {
    return get((int) start);
  }

  @SuppressWarnings("unchecked")
  public T get(int index) {
    return (T) mData[index];
  }

  @Override
  public boolean isEmpty() {
    return false;
  }

  @Override
  public Set<Integer> keySet() {
    Set<Integer> ret = new HashSet<Integer>();

    for (int i : Range.create(mData.length)) {
      if (mData[i] != null) {
        ret.add(i);
      }
    }

    return ret;
  }

  @Override
  public T put(Integer index, T item) {
    return put((int) index, item);
  }

  public T put(int index, T item) {
    mData[index] = item;

    return item;
  }

  // private int getBin(int start) {
  // return start / mBinSize;
  // }

  @Override
  public void putAll(Map<? extends Integer, ? extends T> arg0) {
    // TODO Auto-generated method stub

  }

  @Override
  public T remove(Object arg0) {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public int size() {
    return mData.length;
  }

  @Override
  public Collection<T> values() {
    return null;
  }

  @Override
  public Iterator<Entry<Integer, T>> iterator() {
    return new BinMapIter();
  }

  @Override
  public Entry<Integer, T> first() {
    return iterator().next();
  }
}
