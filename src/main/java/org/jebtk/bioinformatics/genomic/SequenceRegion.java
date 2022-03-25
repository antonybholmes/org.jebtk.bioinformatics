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

import java.util.ArrayList;
import java.util.List;

import org.jebtk.core.text.TextUtils;

/**
 * Describes a region of a genome and its sequence.
 *
 * @author Antony Holmes
 *
 */
public class SequenceRegion extends GenomicRegion {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /** The m sequence. */
  private Sequence mSequence;

  /**
   * Instantiates a new sequence region.
   *
   * @param chr      the chr
   * @param start    the start
   * @param end      the end
   * @param sequence the sequence
   */
  public SequenceRegion(Chromosome chr, int start, int end, Sequence sequence) {
    this(new GenomicRegion(chr, start, end), sequence);
  }

  /**
   * Instantiates a new sequence region.
   *
   * @param region   the region
   * @param sequence the sequence
   */
  public SequenceRegion(GenomicRegion region, String dna) {
    this(region, Sequence.create(dna));
  }

  public SequenceRegion(GenomicRegion region, Sequence sequence) {
    super(region);

    mSequence = sequence;
  }

  /**
   * Gets the region.
   *
   * @return the region
   */
  public Sequence getSequence() {
    return mSequence;
  }

  /**
   * Reverse complement a list of sequences.
   *
   * @param sequences the sequences
   * @return the list
   */
  public static List<SequenceRegion> reverseComplementRegion(List<SequenceRegion> sequences) {
    List<SequenceRegion> ret = new ArrayList<SequenceRegion>(sequences.size());

    for (SequenceRegion sequence : sequences) {
      ret.add(reverseComplement(sequence));
    }

    return ret;
  }

  /**
   * Reverse complement.
   *
   * @param sequence the sequence
   * @return the sequence region
   */
  public static SequenceRegion reverseComplement(SequenceRegion sequence) {
    return new SequenceRegion(GenomicRegion.oppositeStrand(sequence), Sequence.reverseComplement(sequence.mSequence));
  }

  /**
   * Seq to index seq.
   *
   * @param <X>  the generic type
   * @param seqs the seqs
   * @return the char[][]
   */
  public static <X extends SequenceRegion> char[][] seqToIndexSeq(List<X> seqs) {
    char[][] ret = new char[seqs.size()][];

    for (int i = 0; i < seqs.size(); ++i) {
      ret[i] = seqs.get(i).getSequence().toArray();
    }

    return ret;
  }

  /**
   * To index.
   *
   * @param <X>  the generic type
   * @param seqs the seqs
   * @return the byte[][]
   */
  public static <X extends SequenceRegion> byte[][] toIndex(List<X> seqs) {
    byte[][] ret = new byte[seqs.size()][];

    for (int i = 0; i < seqs.size(); ++i) {
      ret[i] = seqs.get(i).getSequence().toIndex();
    }

    return ret;
  }

  /**
   * Return a fasta representation of the sequence
   * 
   * @return
   */
  public String toFasta() {
    return new StringBuilder(">").append(toString()).append(TextUtils.NEW_LINE).append(getSequence()).toString();
  }
}
