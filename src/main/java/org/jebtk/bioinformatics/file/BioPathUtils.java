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
package org.jebtk.bioinformatics.file;

import org.jebtk.core.io.ExtTest;
import org.jebtk.core.io.PathUtils.Ext;

/**
 * The Class BioPathUtils.
 */
public class BioPathUtils {

  /**
   * The Class BioExt.
   */
  public static class BioExt extends Ext {

    /**
     * Bedgraph.
     *
     * @return the ext test
     */
    public ExtTest bedgraph() {
      return type("bedgraph");
    }

    /**
     * Bed.
     *
     * @return the ext test
     */
    public ExtTest bed() {
      return type("bed");
    }

    /**
     * Xlsx.
     *
     * @return the ext test
     */
    public ExtTest xlsx() {
      return type("xlsx");
    }

    /**
     * Xls.
     *
     * @return the ext test
     */
    public ExtTest xls() {
      return type("xls");
    }

    /**
     * Bam.
     *
     * @return the ext test
     */
    public ExtTest bam() {
      return type("bam");
    }

    /**
     * Brt 2.
     *
     * @return the ext test
     */
    public ExtTest brt2() {
      return type("brt2j");
    }

    /**
     * Bvt.
     *
     * @return the ext test
     */
    public ExtTest bvt() {
      return type("bvtj");
    }

    /**
     * Bct.
     *
     * @return the ext test
     */
    public ExtTest bct() {
      return type("bctj");
    }
  }

  /**
   * Ext.
   *
   * @return the bio ext
   */
  public static BioExt ext() {
    return new BioExt();
  }
}
