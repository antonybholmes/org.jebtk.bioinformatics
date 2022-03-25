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
import java.io.IOException;
import java.nio.file.Path;
import java.util.regex.Pattern;

import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.sys.SysUtils;

/**
 * The class Fasta.
 */
public class FastaReader {

  /**
   * The constant FASTA_START.
   */
  private static final String FASTA_START = ">";

  /**
   * The header pattern.
   */
  public static Pattern HEADER_PATTERN = Pattern.compile(">(.+)");

  private BufferedReader mReader = null;

  private String mCurrentName = null;

  private StringBuilder mBuffer = new StringBuilder();

  /**
   * Parses the.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public FastaReader(Path file) throws IOException {
    mReader = FileUtils.newBufferedReader(file);
  }

  public Sequence next() throws IOException {
    String line;

    while ((line = mReader.readLine()) != null) {
      if (Io.isEmptyLine(line)) {
        continue;
      }

      if (line.startsWith(FASTA_START)) {
        Sequence ret = null;

        if (mBuffer.length() > 0) {
          // we were reading a fasta sequence and we have come to the next
          // so output the current

          ret = Sequence.create(mCurrentName, mBuffer.toString());

          // Reset buffer to hold the sequence
          mBuffer.setLength(0);
        }

        // Change the current name otherwise the line read will be wasted
        mCurrentName = line.substring(1);

        // If we found a sequence, return it
        if (ret != null) {
          return ret;
        }
      }

      // Cache the lines
      mBuffer.append(line);
    }

    // Close the reader since we have run out of lines
    // mReader.close();

    if (mBuffer.length() > 0) {
      return Sequence.create(mCurrentName, mBuffer.toString());
    } else {
      return null;
    }
  }

  /**
   * Find the next sequence in the fasta file and put the sequence into the
   * buffer.
   * 
   * @param buffer
   * @return The number of bases read.
   * @throws IOException
   */
  public int next(char[] buffer) throws IOException {
    String line;

    int n = 0;

    char[] buf;

    while ((line = mReader.readLine()) != null) {
      if (Io.isEmptyLine(line)) {
        continue;
      }

      if (line.startsWith(FASTA_START)) {
        // Change the current name otherwise the line read will be wasted
        mCurrentName = line.substring(1);

        // If we found a sequence, return it
        if (n > 0) {
          return n;
        }
      }

      buf = line.toCharArray();

      SysUtils.arraycopy(buf, buffer, n, buf.length);

      n += buf.length;
    }

    // Close the reader since we have run out of lines
    // mReader.close();

    return n;
  }

  /**
   * Returns the current name of the sequence after <code>next()</code> has been
   * called.
   * 
   * @return
   */
  public String currentName() {
    return mCurrentName;
  }

  public void close() throws IOException {
    mReader.close();
  }

  /**
   * Gets the header.
   *
   * @param text the text
   * @return the header
   */
  public static String getHeader(String text) {
    return FASTA_START + text;
  }

}
