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
package org.jebtk.bioinformatics.ext.genepattern;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.math.matrix.Matrix;
import org.jebtk.math.matrix.MatrixParser;

/**
 * The class GctMatrixParser.
 */
public class GctMatrixParser implements MatrixParser {

  /**
   * Sets the. O
   * 
   * @param matrix the matrix
   * @param row    the row
   * @param column the column
   * @param value  the value
   */
  protected void set(Matrix matrix, int row, int column, String value) {

    matrix.set(row, column, Double.parseDouble(value));
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.lib.math.matrix.MatrixParser#parse(java.io.Path)
   */
  @Override
  public DataFrame parse(Path file) throws IOException {
    DataFrame matrix = null;

    BufferedReader reader = FileUtils.newBufferedReader(file);

    String line;

    int r = 0;
    int c = 0;

    List<String> tokens;

    try {
      // skip #1.2

      line = reader.readLine();

      line = reader.readLine();

      tokens = TextUtils.tabSplit(line);

      r = Integer.parseInt(tokens.get(0));
      c = Integer.parseInt(tokens.get(1));

      matrix = GctMatrix.createGctMatrix(r, c);

      List<String> rowNames = new ArrayList<String>();
      List<String> descriptionNames = new ArrayList<String>();

      line = reader.readLine();

      // Look at the columns
      tokens = TextUtils.tabSplit(line);

      String idHeader = tokens.get(0);
      String descriptionHeader = tokens.get(1);

      matrix.setColumnNames(CollectionUtils.subList(tokens, 2));

      int row = 0;

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.tabSplit(TextUtils.removeQuotes(line));

        rowNames.add(tokens.get(0));

        String description = tokens.get(1);

        if (description.length() == 0) {
          description = TextUtils.NA;
        }

        descriptionNames.add(description);

        // the first token is the column name so ignore it
        for (int i = 2; i < tokens.size(); ++i) {
          set(matrix, row, i - 2, tokens.get(i));
        }

        ++row;
      }

      matrix.getIndex().setAnnotation(idHeader, rowNames.toArray());

      matrix.getIndex().setAnnotation(descriptionHeader, descriptionNames.toArray());
    } finally {
      reader.close();
    }

    return matrix;
  }
}
