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

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.math.matrix.DoubleMatrix;

/**
 * Represents a UCSC GCT file in matrix form.
 * 
 * @author Antony Holmes
 *
 */
public class GctMatrix extends DataFrame {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant DESCRIPTION_COLUMN.
   */
  public static final String DESCRIPTION_COLUMN = "Description";

  /**
   * The constant ID_COLUMN.
   */
  private static final String ID_COLUMN = "ID";

  /**
   * The constant VERSION_ID.
   */
  private static final String VERSION_ID = "#1.2";

  /**
   * Instantiates a new gct matrix.
   *
   * @param rows    the rows
   * @param columns the columns
   */
  public GctMatrix(int rows, int columns) {
    super(DoubleMatrix.createDoubleMatrix(rows, columns));
  }

  /**
   * Sets the description names.
   *
   * @param names the new description names
   */
  public void setDescriptionNames(String[] names) {
    getIndex().setAnnotation(DESCRIPTION_COLUMN, names);
  }

  /**
   * Gets the description name.
   *
   * @param i the i
   * @return the description name
   */
  public String getDescriptionName(int i) {
    return getIndex().getText(DESCRIPTION_COLUMN, i);
  }

  /**
   * Gets the description names.
   *
   * @return the description names
   */
  public String[] getDescriptionNames() {
    return getIndex().getText(DESCRIPTION_COLUMN);
  }

  /**
   * Parses the matrix.
   *
   * @param file the file
   * @return the annotation matrix
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static DataFrame parseMatrix(Path file) throws IOException {
    return new GctMatrixParser().parse(file);
  }

  /**
   * Write a simple expression matrix in GCT format.
   *
   * @param <T>    the generic type
   * @param matrix the matrix
   * @param file   the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static <T> void writeGctMatrix(DataFrame matrix, Path file) throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    try {
      writer.write(VERSION_ID);
      writer.newLine();
      writer.write(Integer.toString(matrix.getRows()));
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(Integer.toString(matrix.getCols()));
      writer.newLine();

      writer.write(ID_COLUMN);
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(DESCRIPTION_COLUMN);

      for (int i = 0; i < matrix.getCols(); ++i) {
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(matrix.getColumnName(i));
      }

      writer.newLine();

      List<String> names = matrix.getIndex().getNames();

      for (int i = 0; i < matrix.getRows(); ++i) {
        writer.write(matrix.getIndex().getText(names.get(0), i));
        writer.write(TextUtils.TAB_DELIMITER);

        writer.write(matrix.getIndex().getText(names.get(names.size() - 1), i));

        for (int j = 0; j < matrix.getCols(); ++j) {
          writer.write(TextUtils.TAB_DELIMITER);

          writer.write(Double.toString(matrix.getValue(i, j)));
        }

        writer.newLine();
      }
    } finally {
      writer.close();
    }
  }

  /**
   * Creates the gct matrix.
   *
   * @param rows the rows
   * @param cols the cols
   * @return the annotation matrix
   */
  public static DataFrame createGctMatrix(int rows, int cols) {
    return new GctMatrix(rows, cols);
  }
}
