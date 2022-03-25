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

import com.fasterxml.jackson.annotation.JsonGetter;

/**
 * The Class Species.
 */
public class Species extends Type {

  /** The m scientific name. */
  private String mScientificName;

  /**
   * Instantiates a new species.
   *
   * @param name           the name
   * @param scientificName the scientific name
   */
  public Species(String name, String scientificName) {
    this(-1, name, scientificName);
  }

  /**
   * Instantiates a new species.
   *
   * @param name the name
   */
  public Species(String name) {
    this(-1, name);
  }

  /**
   * Instantiates a new species.
   *
   * @param id   the id
   * @param name the name
   */
  public Species(int id, String name) {
    this(id, name, name);
  }

  /**
   * Instantiates a new species.
   *
   * @param id             the id
   * @param name           the name
   * @param scientificName the scientific name
   */
  public Species(int id, String name, String scientificName) {
    super(id, name);

    mScientificName = scientificName;
  }

  /**
   * Gets the scientific name.
   *
   * @return the scientific name
   */
  @JsonGetter("sn")
  public String getScientificName() {
    return mScientificName;
  }
}
