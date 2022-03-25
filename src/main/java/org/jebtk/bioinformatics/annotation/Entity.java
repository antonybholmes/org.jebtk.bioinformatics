/**
 * Copyright 2017 Antony Holmes
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
package org.jebtk.bioinformatics.annotation;

import org.jebtk.core.IdProperty;

import com.fasterxml.jackson.annotation.JsonGetter;

/**
 * The Class Entity.
 */
public class Entity implements IdProperty, Comparable<Entity> {

  /** The m id. */
  public final int mId;

  /**
   * Instantiates a new entity.
   *
   * @param id the id
   */
  public Entity(int id) {
    mId = id;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.common.IdProperty#getId()
   */
  @Override
  @JsonGetter("id")
  public int getId() {
    return mId;
  }

  @Override
  public int compareTo(Entity t) {
    if (mId > t.mId) {
      return 1;
    } else if (mId < t.mId) {
      return -1;
    } else {
      return 0;
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof Entity) {
      return compareTo((Entity) o) == 0;
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
    return mId;
  }
}
