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
package org.jebtk.bioinformatics.dna;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.core.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Fast search of genome sequence Paths to get get actual genomic data.
 *
 * @author Antony Holmes
 */
public class GZipSequenceReader extends ChrSequenceReader {

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(GZipSequenceReader.class);

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public GZipSequenceReader(Path directory) {
    super(directory);
  }

  @Override
  public String getName() {
    return "txt-gz";
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.genome.GenomeAssembly#getSequence(edu.
   * columbia.rdf.lib.bioinformatics.genome.GenomicRegion, boolean,
   * edu.columbia.rdf.lib.bioinformatics.genome.RepeatMaskType)
   */
  @Override
  public final SequenceRegion getSequence(Genome genome, GenomicRegion region, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {
    Chromosome chr = region.getChr();

    if (!mFileMap.containsKey(region.getChr())) {
      mFileMap.put(chr, mFile.resolve(chr + ".txt.gz"));
    }

    return getSequence(mFileMap.get(chr), region);
  }

  /**
   * Returns the sequence from a region of a chromosome.
   *
   * @param file   the file
   * @param region the region
   * @return the sequence
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private final SequenceRegion getSequence(Path file, GenomicRegion region) throws IOException {
    LOG.debug("Extract sequence for {} from {}...", region.toString(), file);

    int s = region.getStart() - 1;
    int e = region.getEnd() - 1;

    // PathInputStream in = new PathInputStream(Path);
    InputStream in = FileUtils.newGzipInputStream(file); // new
                                                         // GZIPInputStream(new
                                                         // PathInputStream(Path),
                                                         // 65536);

    int l = e - s + 1;

    byte[] cbuf = new byte[l];

    SequenceRegion sequence = null;

    try {
      in.skip(s);

      int bytesRead = in.read(cbuf);

      if (bytesRead == l) {
        sequence = new SequenceRegion(region, Sequence.create(new String(cbuf)));
      } else {
        System.err.println(file + " " + region + " " + bytesRead + " " + l + " " + new String(cbuf));
        System.exit(0);
      }
    } finally {
      in.close();
    }

    return sequence;
  }
}
