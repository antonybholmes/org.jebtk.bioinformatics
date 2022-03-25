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
import org.jebtk.math.matrix.Matrix;
import org.jebtk.math.matrix.MixedMatrix;

/**
 * The class ResMatrix.
 */
public class ResMatrix extends DataFrame {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant DESCRIPTION_COLUMN.
   */
  public static final String DESCRIPTION_COLUMN = "Description";

  /**
   * The constant ACCESSION_COLUMN.
   */
  public static final String ACCESSION_COLUMN = "Accession";

  /**
   * The constant SAMPLE_DESCRIPTIONS.
   */
  private static final String SAMPLE_DESCRIPTIONS = "Sample Description";

  /**
   * Instantiates a new res matrix.
   *
   * @param rows    the rows
   * @param columns the columns
   */
  public ResMatrix(int rows, int columns) {
    this(new MixedMatrix(rows, columns));
  }

  /**
   * Instantiates a new res matrix.
   *
   * @param matrix the matrix
   */
  private ResMatrix(Matrix matrix) {
    super(matrix);
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
   * Sets the accession names.
   *
   * @param names the new accession names
   */
  public void setAccessionNames(String[] names) {
    getIndex().setAnnotation(ACCESSION_COLUMN, names);
  }

  /**
   * Gets the accession name.
   *
   * @param i the i
   * @return the accession name
   */
  public String getAccessionName(int i) {
    return getIndex().getText(ACCESSION_COLUMN, i);
  }

  /**
   * Gets the accession names.
   *
   * @return the accession names
   */
  public String[] getAccessionNames() {
    return getIndex().getText(ACCESSION_COLUMN);
  }

  /**
   * Sets the sample descriptions.
   *
   * @param names the new sample descriptions
   */
  public void setSampleDescriptions(String[] names) {
    getColumnHeader().setAnnotation(SAMPLE_DESCRIPTIONS, names);
  }

  /**
   * Gets the sample description.
   *
   * @param i the i
   * @return the sample description
   */
  public String getSampleDescription(int i) {
    return getIndex().getText(SAMPLE_DESCRIPTIONS, i);
  }

  /**
   * Gets the sample descriptions.
   *
   * @return the sample descriptions
   */
  public String[] getSampleDescriptions() {
    return getIndex().getText(SAMPLE_DESCRIPTIONS);
  }

  /**
   * Parses the matrix.
   *
   * @param file the file
   * @return the annotation matrix
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static DataFrame parseMatrix(Path file) throws IOException {
    return parseMatrix(file, false);
  }

  /**
   * Parses the matrix.
   *
   * @param file         the file
   * @param keepCallCols the keep call cols
   * @return the annotation matrix
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static DataFrame parseMatrix(Path file, boolean keepCallCols) throws IOException {
    return new ResMatrixParser(keepCallCols).parse(file);
  }

  /**
   * Write a simple expression matrix in GCT format.
   *
   * @param <T>    the generic type
   * @param matrix the matrix
   * @param file   the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static <T> void writeResMatrix(DataFrame matrix, Path file) throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    try {
      writer.write(DESCRIPTION_COLUMN);
      writer.write(TextUtils.TAB_DELIMITER);
      writer.write(ACCESSION_COLUMN);

      for (int i = 0; i < matrix.getCols(); ++i) {
        if (i > 0) {
          writer.write(TextUtils.TAB_DELIMITER);
        }

        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(matrix.getColumnName(i));
      }

      writer.newLine();

      // Insert a blank line

      // for (int i = 0; i < matrix.getColumnCount(); ++i) {
      // writer.write(TextUtils.TAB_DELIMITER);
      // writer.write(matrix.getSampleDescription(i));
      // }

      writer.newLine();

      // Write row count
      writer.write(Integer.toString(matrix.getRows()));
      writer.newLine();

      List<String> names = matrix.getIndex().getNames();

      for (int i = 0; i < matrix.getRows(); ++i) {
        writer.write(matrix.getIndex().getText(names.get(0), i));
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(matrix.getIndex().getText(names.get(names.size() - 1), i));

        for (int j = 0; j < matrix.getCols(); ++j) {
          if (i > 0) {
            writer.write(TextUtils.TAB_DELIMITER);
          }

          writer.write(TextUtils.TAB_DELIMITER);
          writer.write(formatTextValue(matrix.getText(i, j)));
        }

        writer.newLine();
      }
    } finally {
      writer.close();
    }
  }

  /**
   * Format text value.
   *
   * @param <T>   the generic type
   * @param value the value
   * @return the string
   */
  public static <T> String formatTextValue(T value) {
    if (value == null) {
      return TextUtils.EMPTY_STRING;
    }

    return value.toString();
  }

  /**
   * Creates the mixed res matrix.
   *
   * @param rows    the rows
   * @param columns the columns
   * @return the res matrix
   */
  public static ResMatrix createMixedResMatrix(int rows, int columns) {
    return new ResMatrix(MixedMatrix.createMixedMatrix(rows, columns));
  }

  /**
   * Creates the res matrix.
   *
   * @param rows    the rows
   * @param columns the columns
   * @return the res matrix
   */
  public static ResMatrix createResMatrix(int rows, int columns) {
    return new ResMatrix(DoubleMatrix.createDoubleMatrix(rows, columns));
  }
}
