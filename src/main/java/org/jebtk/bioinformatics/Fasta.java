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
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * The class Fasta.
 */
public class Fasta {

  /**
   * The constant FASTA_START.
   */
  private static final String FASTA_START = ">";

  /**
   * The header pattern.
   */
  public static Pattern HEADER_PATTERN = Pattern.compile(">(.+)");

  /**
   * Parses the.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static final List<Sequence> parse(Path file) throws IOException {
    List<Sequence> sequences = new ArrayList<Sequence>();

    // System.out.println(file.toString());

    BufferedReader reader = FileUtils.newBufferedReader(file);

    String line;

    Matcher fastaHeaderMatcher;

    String name = null;
    StringBuilder buffer = null;

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        fastaHeaderMatcher = HEADER_PATTERN.matcher(line);

        if (fastaHeaderMatcher.find()) {
          if (buffer != null) {
            sequences.add(Sequence.create(name, buffer.toString()));
          }

          name = fastaHeaderMatcher.group(1);

          buffer = new StringBuilder();
        } else {
          buffer.append(line);
        }
      }
    } finally {
      reader.close();
    }

    // Add the last sequence read
    sequences.add(Sequence.create(name, buffer.toString()));

    return sequences;
  }

  /**
   * Write.
   *
   * @param file     the file
   * @param sequence the sequence
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(Path file, Sequence sequence) throws IOException {
    write(file, CollectionUtils.asList(sequence));
  }

  /**
   * Write a series of fasta sequences to a file.
   *
   * @param file      the file
   * @param sequences the sequences
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(Path file, Collection<Sequence> sequences) throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    try {
      for (Sequence sequence : sequences) {
        writer.write(getHeader(sequence.getName()));
        writer.newLine();
        writer.write(sequence.toString());
        writer.newLine();
      }
    } finally {
      writer.close();
    }
  }

  /**
   * Write.
   *
   * @param file     the file
   * @param sequence the sequence
   * @param width    the width
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(Path file, Sequence sequence, int width) throws IOException {
    write(file, CollectionUtils.asList(sequence), width);
  }

  /**
   * Write.
   *
   * @param file      the file
   * @param sequences the sequences
   * @param width     the width
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(Path file, Collection<Sequence> sequences, int width) throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    try {
      for (Sequence sequence : sequences) {
        writer.write(getHeader(sequence.getName()));
        writer.newLine();

        List<String> lines = TextUtils.breakApart(sequence.toString(), width);

        for (String line : lines) {
          writer.write(line);
          writer.newLine();
        }
      }
    } finally {
      writer.close();
    }
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
