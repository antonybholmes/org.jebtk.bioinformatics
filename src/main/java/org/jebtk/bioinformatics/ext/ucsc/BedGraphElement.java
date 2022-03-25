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

import java.io.IOException;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.ChromosomeService;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.core.text.TextUtils;

/**
 * Genomic region plus value.
 * 
 * @author Antony Holmes
 *
 */
public class BedGraphElement extends GenomicElement {

  private static final long serialVersionUID = 1L;

  /**
   * Create a new region. Coordinates are one based.
   *
   * @param chromosome the chromosome
   * @param start      the start
   * @param end        the end
   * @param value      the value
   */
  public BedGraphElement(GenomicType type, GenomicRegion r, double value) {
    super(type, r);

    setProperty("value", value);
  }

  /**
   * Gets the value.
   *
   * @return the value
   */
  public double getValue() {
    return getDouble("value");
  }

  /*
   * (non-Javadoc)
   * 
   * @see edu.columbia.rdf.lib.bioinformatics.external.ucsc.UCSCTrackRegion#
   * formattedTxt(java.lang.Appendable)
   */
  @Override
  public void formattedTxt(Appendable buffer) throws IOException {
    super.formattedTxt(buffer);
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Double.toString(getValue()));
  }

  public static BedGraphElement parse(Genome genome, String line) {
    return parse(GenomicType.REGION, genome, line);
  }

  /**
   * Parses the.
   *
   * @param line the line
   * @return the bed graph region
   */
  public static BedGraphElement parse(GenomicType type, Genome genome, String line) {
    List<String> tokens = TextUtils.fastSplit(line, TextUtils.TAB_DELIMITER);

    // convert first part to chromosome (replacing x,y and m) {
    Chromosome chr = ChromosomeService.getInstance().chr(genome, tokens.get(0));

    // Per UCSC convention, coordinates are zero based in the file
    // but we shall convert them to one based for consistency
    int start = Integer.parseInt(tokens.get(1)) + 1;
    int end = Integer.parseInt(tokens.get(2));
    double value = Double.parseDouble(tokens.get(3));

    return new BedGraphElement(type, new GenomicRegion(chr, start, end), value);
  }
}
