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
package org.jebtk.bioinformatics.ext.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.ProbeGene;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

/**
 * The class ChipFile.
 */
public class ChipFile {

  /**
   * Instantiates a new chip file.
   */
  private ChipFile() {
    // do nothing
  }

  /**
   * Parses the chip file.
   *
   * @param file the file
   * @return the map
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Map<String, ProbeGene> parseChipFile(Path file) throws IOException {
    System.err.println("Parsing chip file " + file);

    Map<String, ProbeGene> probeGeneMap = new HashMap<String, ProbeGene>();

    BufferedReader reader = FileUtils.newBufferedReader(file);

    String line;

    Splitter split = Splitter.onTab();

    try {
      // skip header
      List<String> header = split.text(reader.readLine());

      int idCol = TextUtils.findFirst(header, "Probe Set ID");
      int symbolCol = TextUtils.findFirst(header, "Gene Symbol");
      int titleCol = TextUtils.findFirst(header, "Gene Title");

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        List<String> tokens = split.text(line);

        String probe = tokens.get(idCol);

        ProbeGene gene;

        if (titleCol != -1) {
          gene = new ProbeGene(tokens.get(idCol), tokens.get(symbolCol), tokens.get(titleCol));
        } else {
          gene = new ProbeGene(tokens.get(idCol), tokens.get(symbolCol));
        }

        probeGeneMap.put(probe, gene);

        // System.err.println(probe + " " + gene);
      }

    } finally {
      reader.close();
    }

    return probeGeneMap;
  }
}
