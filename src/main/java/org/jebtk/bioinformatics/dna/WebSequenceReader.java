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
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.DataSource;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomeService;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceReader;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.core.http.URLPath;
import org.jebtk.core.json.Json;
import org.jebtk.core.json.JsonParser;

// TODO: Auto-generated Javadoc
/**
 * Maintains a connection to a caArray server.
 *
 * @author Antony Holmes
 */
public class WebSequenceReader extends SequenceReader {
  /**
   * The member url.
   */
  private URLPath mUrl;

  /**
   * The member parser.
   */
  private JsonParser mParser;

  /** The m genomes url. */
  private URLPath mGenomesUrl;

  /** The m dna url. */
  private URLPath mDnaUrl;

  /**
   * Instantiates a new genome assembly web.
   *
   * @param url the url
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public WebSequenceReader(URL url) throws IOException {
    mUrl = URLPath.fromUrl(url);

    mDnaUrl = mUrl.join("seq");
    mGenomesUrl = mUrl.join("genomes");

    mParser = new JsonParser();
  }

  @Override
  public String getName() {
    return "http";
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
  public SequenceRegion getSequence(Genome genome, GenomicRegion region, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {
    URL url;

    try {
      URLPath tmpUrl = mDnaUrl.param("n", genome.getName().toLowerCase()).param("a", genome.getAssembly())
          // .param("t", region.getGenome().getTrack())
          .param("chr", region.getChr().toString()).param("s", region.getStart()).param("e", region.getEnd())
          .param("strand", "s").param("lc", displayUpper ? "0" : "1");

      switch (repeatMaskType) {
      case UPPERCASE:
        tmpUrl = tmpUrl.param("mask", "u");
        break;
      case N:
        tmpUrl = tmpUrl.param("mask", "n");
        break;
      default:
        tmpUrl = tmpUrl.param("mask", "l");
        break;
      }

      url = tmpUrl.toURL();

      System.err.println(url);

      Json json = mParser.parse(url);

      String dna = json.get("seq").getString();

      SequenceRegion ret = new SequenceRegion(region, Sequence.create(dna));

      return ret;
    } catch (MalformedURLException e) {
      e.printStackTrace();
    }

    return null;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genome.GenomeAssembly#getGenomes()
   */
  @Override
  public List<Genome> getGenomes() throws IOException {

    List<Genome> ret = new ArrayList<Genome>(100);

    URL url;

    try {
      url = mGenomesUrl.toURL();

      System.err.println(url);

      Json json = mParser.parse(url);
      Json e;

      for (int i = 0; i < json.size(); ++i) {
        e = json.get(i);

        System.err.println("web " + e.getString());

        ret.add(GenomeService.getInstance().get(e.getString("name"), e.getString("assembly")));

        // GenomeService.getInstance().guessGenome(json.get(i).getString()));
      }
    } catch (MalformedURLException e) {
      e.printStackTrace();
    }

    return ret;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genome.GenomeAssembly#getDataSource()
   */
  @Override
  public DataSource getDataSource() {
    return DataSource.WEB;
  }
}
