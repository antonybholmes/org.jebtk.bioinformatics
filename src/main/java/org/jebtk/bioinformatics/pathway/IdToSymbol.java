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
package org.jebtk.bioinformatics.pathway;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * Convert between gene ids/symbols.
 */
public class IdToSymbol {

  /**
   * The member ref seq map.
   */
  private Map<String, String> mId1Map = new HashMap<String, String>();

  /**
   * The member symbol map.
   */
  private Set<String> mSymbolSet = new TreeSet<String>();

  /**
   * Create a conversion tool. Gene symbols
   *
   * @param reader the reader
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public IdToSymbol(BufferedReader reader) throws IOException {

    try {
      // Skip header
      reader.readLine();

      String line;

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        List<String> tokens = TextUtils.tabSplit(line);

        String id1 = tokens.get(0);
        String id2 = tokens.get(1);
        String symbol = tokens.get(2);

        mId1Map.put(id1.toUpperCase(), symbol);
        mId1Map.put(id2.toUpperCase(), symbol);
        mId1Map.put(symbol.toUpperCase(), symbol);

        mSymbolSet.add(symbol);
      }
    } finally {
      reader.close();
    }
  }

  /**
   * Convert.
   *
   * @param id the id
   * @return the string
   */
  public String convert(String id) {
    return mId1Map.get(id.toUpperCase());
  }

  /**
   * Gets the symbol count.
   *
   * @return the symbol count
   */
  public int getSymbolCount() {
    return mSymbolSet.size();
  }
}
