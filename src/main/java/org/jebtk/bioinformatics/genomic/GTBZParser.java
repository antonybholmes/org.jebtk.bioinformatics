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
package org.jebtk.bioinformatics.genomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jebtk.core.io.FileUtils;

/**
 * Genes lookup to m.
 *
 * @author Antony Holmes
 */
public class GTBZParser extends GTB2Parser {

  public GTBZParser() {
    // Do nothing
  }

  public GTBZParser(GeneParser parser) {
    super(parser);
  }

  @Override
  public void parse(Path file, Genome genome, GenesDB genes) throws IOException {

    final ZipFile zipFile = new ZipFile(file.toFile());

    try {
      final Enumeration<? extends ZipEntry> entries = zipFile.entries();

      while (entries.hasMoreElements()) {
        final ZipEntry entry = entries.nextElement();

        if (entry.getName().contains("gtb2")) {
          BufferedReader reader = FileUtils.newBufferedReader(zipFile, entry);

          try {
            parse(file, reader, genome, genes);
          } finally {
            reader.close();
          }
        }
      }
    } finally {
      zipFile.close();
    }
  }

  @Override
  public void parse(Path file, Genome genome, Chromosome chr, GenesDB genes) throws IOException {

    final ZipFile zipFile = new ZipFile(file.toFile());

    try {
      ZipEntry entry = zipFile.getEntry(chr.getName() + ".gtb2");

      if (entry != null) {
        BufferedReader reader = FileUtils.newBufferedReader(zipFile, entry);

        try {
          parse(file, reader, genome, genes);
        } finally {
          reader.close();
        }
      }
    } finally {
      zipFile.close();
    }
  }

  @Override
  public GeneParser create(GeneParser parser) {
    GTBZParser ret = new GTBZParser(parser);

    return ret;
  }
}