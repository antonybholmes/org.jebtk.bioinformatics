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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.core.Resources;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Server for genome feature annotations.
 *
 * @author Antony Holmes
 *
 */
public class CytobandsService {

  /**
   * The constant INSTANCE.
   */
  private static final CytobandsService INSTANCE = new CytobandsService();

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(CytobandsService.class);

  /**
   * The constant DEFAULT_RES.
   */
  public static final String DEFAULT_HG18_RES = "res/ucsc_cytobands_hg18.txt.gz";

  /**
   * The constant DEFAULT_FILE.
   */
  public static final File DEFAULT_HG18_FILE = new File(DEFAULT_HG18_RES);

  /**
   * The constant DEFAULT_RES.
   */
  public static final String DEFAULT_HG19_RES = "res/ucsc_cytobands_hg19.txt.gz";

  /**
   * The constant DEFAULT_FILE.
   */
  public static final File DEFAULT_HG19_FILE = new File(DEFAULT_HG19_RES);

  /**
   * Gets the single instance of Cytobands.
   *
   * @return single instance of Cytobands
   */
  public static final CytobandsService getInstance() {
    return INSTANCE;
  }

  /**
   * The member cytobands map.
   */
  private Map<Genome, Cytobands> mCytobandsMap = new HashMap<Genome, Cytobands>();

  /**
   * Instantiates a new cytobands.
   */
  public CytobandsService() {
    // do nothing
  }

  /**
   * Load internal hg 18.
   *
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public final void loadInternalHg18() throws IOException {
    LOG.info("Parsing {}...", DEFAULT_HG18_RES);

    load(Genome.HG18, Resources.getResGzipReader(DEFAULT_HG18_RES));
  }

  /**
   * Load the default lib from an external file.
   *
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public final void loadRefSeqHg18() throws IOException {
    LOG.info("Parsing {}...", DEFAULT_HG18_FILE);

    load(Genome.HG18, Resources.getGzipReader(DEFAULT_HG18_FILE));
  }

  /**
   * Load the default lib from the jar resource.
   *
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public final void loadInternalHg19() throws IOException {
    LOG.info("Parsing {}...", DEFAULT_HG19_RES);

    load(Genome.HG19, Resources.getResGzipReader(DEFAULT_HG19_RES));
  }

  /**
   * Load the default lib from an external file.
   *
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public final void loadRefSeqHg19() throws IOException {
    LOG.info("Parsing {}...", DEFAULT_HG19_FILE);

    load(Genome.HG19, Resources.getGzipReader(DEFAULT_HG19_FILE));
  }

  /**
   * Load.
   *
   * @param genome the genome
   * @param reader the reader
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public void load(Genome genome, BufferedReader reader) throws IOException {
    load(genome, new Cytobands(genome, reader));
  }

  /**
   * Load.
   *
   * @param genome    the genome
   * @param cytobands the cytobands
   */
  public void load(Genome genome, Cytobands cytobands) {
    mCytobandsMap.put(genome, cytobands);
  }

  /**
   * Gets the cytobands.
   *
   * @param genome the genome
   * @return the cytobands
   */
  public Cytobands getCytobands(Genome genome) {
    System.err.println("cyoservicekeys:" + genome);
    System.err.println("cyoservicekeys:" + mCytobandsMap.keySet());
    
    System.err.println("fgfdg"+ mCytobandsMap.get(genome));
    
    return mCytobandsMap.get(genome);
  }

}