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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * Allows multiple logs to be agglomerated so a message can be fire to multiple
 * logs.
 *
 * @author Antony Holmes
 *
 */
public class AminoAcids implements Iterable<AminoAcid> {

  /**
   * The constant INSTANCE.
   */
  private static final AminoAcids INSTANCE = new AminoAcids();

  /**
   * Gets the single instance of AminoAcids.
   *
   * @return single instance of AminoAcids
   */
  public static final AminoAcids getInstance() {
    return INSTANCE;
  }

  /**
   * The amino acids.
   */
  private Map<String, AminoAcid> aminoAcids = new HashMap<String, AminoAcid>();

  /**
   * The codon acid map.
   */
  private Map<Codon, AminoAcid> codonAcidMap = new HashMap<Codon, AminoAcid>();

  /**
   * Instantiates a new amino acids.
   */
  private AminoAcids() {
    // do nothing
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<AminoAcid> iterator() {
    return aminoAcids.values().iterator();
  }

  /**
   * Gets the amino acid by name.
   *
   * @param name the name
   * @return the amino acid by name
   */
  public AminoAcid getAminoAcidByName(String name) {
    return aminoAcids.get(name);
  }

  /**
   * Gets the amino acid by codon.
   *
   * @param codon the codon
   * @return the amino acid by codon
   */
  public AminoAcid getAminoAcidByCodon(String codon) {
    return getAminoAcidByCodon(new Codon(codon));
  }

  /**
   * Gets the amino acid by codon.
   *
   * @param codon the codon
   * @return the amino acid by codon
   */
  public AminoAcid getAminoAcidByCodon(Codon codon) {
    return codonAcidMap.get(codon);
  }

  /**
   * Loads a tab delimited list of amino acids. Line 1 is a header Column 1 Amino
   * Acid name Column 2 Short Code Column 3 codons, semi-colon separated.
   *
   * @param file the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public void load(File file) throws IOException {
    aminoAcids.clear();

    BufferedReader reader = new BufferedReader(new FileReader(file));

    String line;
    List<String> tokens;

    try {
      // skip header
      reader.readLine();

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.fastSplit(line, TextUtils.TAB_DELIMITER);

        Codon codon = new Codon(tokens.get(0));
        String name = tokens.get(1);
        String abbreviation = tokens.get(2);
        char letter = tokens.get(3).charAt(0);

        // If the amino acid is alread in the db, update
        // else create a new entry

        if (!aminoAcids.containsKey(name)) {
          aminoAcids.put(name, new AminoAcid(name, abbreviation, letter));
        }

        aminoAcids.get(name).addCodon(codon);

        codonAcidMap.put(codon, aminoAcids.get(name));
      }
    } finally {
      reader.close();
    }
  }
}
