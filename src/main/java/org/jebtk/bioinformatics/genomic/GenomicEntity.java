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

import org.jebtk.core.NameGetter;

import com.fasterxml.jackson.annotation.JsonIgnore;

// TODO: Auto-generated Javadoc
/**
 * Represents a genomic entity such as a gene or exon. An entity is a genomic
 * region with a collection of annotations and sub entities to describe a
 * genomic feature such as a gene, transcript or exon. It is designed for use
 * with different annotation databases such as Refseq and Ensembl where
 * different ids and nomenclature are used.
 */
public class GenomicEntity extends GenomicElement implements NameGetter {
  private static final long serialVersionUID = 1L;

  /** The Constant SYMBOL_TYPE. */
  // public static final String SYMBOL_TYPE = "symbol";

  public static final String GENE_ID = "gene_id";

  /** The Constant GENE_NAME_TYPE. */
  public static final String GENE_NAME = "gene_name";

  /** The Constant TRANSCRIPT_ID_TYPE. */
  public static final String TRANSCRIPT_ID = "transcript_id";

  /** The Constant REFSEQ_TYPE. */
  public static final String REFSEQ_ID = "refseq_id";

  /** The Constant ENTREZ_TYPE. */
  public static final String ENTREZ_ID = "entrez_id";

  /** The Constant ENSEMBL_TYPE. */
  public static final String ENSEMBL_ID = "ensembl_id";

  public static final String UTR_5P = "utr_5p";

  public static final String UTR_3P = "utr_3p";

  /**
   * Instantiates a new genomic entity.
   *
   * @param type the type
   * @param l    the l
   */
  public GenomicEntity(GenomicType type, GenomicRegion l) {
    super(type, l);
  }

  /**
   * Instantiates a new genomic entity.
   *
   * @param type the type
   * @param l    the l
   * @param s    the s
   */
  public GenomicEntity(GenomicType type, GenomicRegion l, Strand s) {
    super(type, l, s);
  }

  public GenomicEntity(GenomicType type, Genome genome, Chromosome chr, int start, int end) {
    this(type, genome, chr, start, end, Strand.SENSE);
  }

  public GenomicEntity(GenomicType type, Genome genome, Chromosome chr, int start, int end, Strand strand) {
    super(type, chr, start, end, strand);
  }

  /**
   * Set a named property of the gene, such as a gene symbol or entrez id.
   *
   * @param name  the type
   * @param value the name
   * @return the genomic entity
   */
  public GenomicElement setProperty(String name, String value) {

    String lc = name.toLowerCase();

    if (lc.contains("symbol")) {
      super.setProperty(GENE_NAME, value);
    }

    if (lc.contains("transcript")) {
      return super.setProperty(TRANSCRIPT_ID, value);
    }

    return super.setProperty(lc, value);
  }

  /**
   * Gets the ref seq.
   *
   * @return the ref seq
   */
  @JsonIgnore
  public String getRefSeq() {
    return getProperty(REFSEQ_ID).toString();
  }

  /**
   * Gets the entrez.
   *
   * @return the entrez
   */
  @JsonIgnore
  public String getEntrez() {
    return getProperty(ENTREZ_ID);
  }

  /**
   * Gets the gene symbol. Equivalent to {@code getGeneName()}.
   *
   * @return the gene symbol.
   */
  @JsonIgnore
  public String getSymbol() {
    return getGeneName();
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.core.NameProperty#getName()
   */
  @JsonIgnore
  @Override
  public String getName() {
    return getGeneName();
  }

  /**
   * Gets the gene id. Equivalent to {@code getGeneId()}.
   *
   * @return the id
   */
  @JsonIgnore
  public String getId() {
    return getGeneId();
  }

  /**
   * Gets the gene id.
   *
   * @return the gene id
   */
  @JsonIgnore
  public String getGeneId() {
    return getProperty(GENE_ID);
  }

  /**
   * Gets the gene name.
   *
   * @return the gene name
   */
  @JsonIgnore
  public String getGeneName() {
    return getProperty(GENE_NAME);
  }

  /**
   * Sets the symbol.
   *
   * @param name the name
   * @return the genomic entity
   */
  public GenomicEntity setSymbol(String name) {
    return setGeneName(name);
  }

  /**
   * Sets the gene name.
   *
   * @param name the name
   * @return the genomic entity
   */
  public GenomicEntity setGeneName(String name) {
    return (GenomicEntity) setProperty(GENE_NAME, name);
  }

  /**
   * Sets the refseq.
   *
   * @param name the name
   * @return the gene
   */
  public GenomicEntity setRefseq(String name) {
    return (GenomicEntity) setProperty(REFSEQ_ID, name);
  }

  /**
   * Sets the entrez.
   *
   * @param name the name
   * @return the gene
   */
  public GenomicEntity setEntrez(String name) {
    return (GenomicEntity) setProperty(ENTREZ_ID, name);
  }

  /**
   * Gets the transcript id.
   *
   * @return the transcript id
   */
  @JsonIgnore
  public String getTranscriptId() {
    return getProperty(TRANSCRIPT_ID).toString();
  }

  /**
   * Sets the transcript id.
   *
   * @param name the name
   * @return the genomic entity
   */
  public GenomicEntity setTranscriptId(String name) {
    return (GenomicEntity) setProperty(TRANSCRIPT_ID, name);
  }

}
