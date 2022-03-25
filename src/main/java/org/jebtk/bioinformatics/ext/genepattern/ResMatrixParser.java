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
import java.util.List;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.stream.Stream;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.math.matrix.Matrix;
import org.jebtk.math.matrix.MatrixParser;

/**
 * Parses a simple txt matrix.
 *
 * @author Antony Holmes
 */
public class ResMatrixParser implements MatrixParser {

  /** The m keep call col. */
  private boolean mKeepCallCol;

  /**
   * Instantiates a new res matrix parser.
   *
   * @param keepCallCol the keep call col
   */
  public ResMatrixParser(boolean keepCallCol) {
    mKeepCallCol = keepCallCol;
  }

  /**
   * Sets the.
   *
   * @param matrix      the matrix
   * @param row         the row
   * @param column      the column
   * @param value       the value
   * @param keepCallCol the keep call col
   */
  private static void set(Matrix matrix, int row, int column, String value, boolean keepCallCol) {
    if (keepCallCol) {
      if (column % 2 == 0) {
        matrix.set(row, column, Double.parseDouble(value));
      } else {
        // The second column in each group is the call column so
        // no need to parse
        matrix.set(row, column, value);
      }
    } else {
      matrix.set(row, column, Double.parseDouble(value));
    }
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

    List<String> tokens;

    int rows = 0;
    int columns = -1;

    try {
      // skip first line
      // line = reader.readLine();

      line = reader.readLine();
      tokens = TextUtils.tabSplit(line);

      // columns is cols - 2 annotation columns then divided
      // by two to account for present/absent
      columns = tokens.size() - 2;

      if (!mKeepCallCol) {
        // If we are not keeping the call column, we only need half
        // as many columns
        columns /= 2;
      }

      // skip second header
      reader.readLine();

      String s = reader.readLine();

      if (TextUtils.isInt(s)) {
        rows = Integer.parseInt(s);
      } else {
        // Since we read the line to test if it was the number of
        // rows, we start at 1 for the count since the line was in
        // fact a valid row
        rows = 1;

        while ((line = reader.readLine()) != null) {
          if (Io.isEmptyLine(line)) {
            continue;
          }

          ++rows;
        }
      }
    } finally {
      reader.close();
    }

    if (mKeepCallCol) {
      matrix = ResMatrix.createMixedResMatrix(rows, columns);
    } else {
      matrix = ResMatrix.createResMatrix(rows, columns);
    }

    reader = FileUtils.newBufferedReader(file);

    try {
      // skip header
      // reader.readLine();

      // add column names
      line = reader.readLine();

      Stream<String> stream = Stream.of(TextUtils.tabSplit(line)).skip(2).jump(2);

      if (mKeepCallCol) {
        stream = stream.replicate(2).asString().append(" call", 1, 2);

        // tokens = CollectionUtils.replicate(tokens, 2);
        // tokens = TextUtils.append(tokens, " call", 2);
      }

      tokens = stream.toList();

      matrix.setColumnNames(tokens);

      // skip second header
      reader.readLine();

      int row = 0;

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        tokens = TextUtils.tabSplit(line);

        matrix.setRowName(row, tokens.get(1));
        matrix.getIndex().setAnnotation("Description", row, tokens.get(0));

        stream = Stream.of(tokens).skip(2);

        if (!mKeepCallCol) {
          stream = stream.jump(2);
        }

        tokens = stream.toList();

        // the first token is the column name so ignore it
        for (int i = 0; i < tokens.size(); ++i) {
          set(matrix, row, i, tokens.get(i), mKeepCallCol);
        }

        ++row;
      }
    } finally {
      reader.close();
    }

    return matrix;
  }
}
