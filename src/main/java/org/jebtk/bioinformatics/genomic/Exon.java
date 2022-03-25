package org.jebtk.bioinformatics.genomic;

public class Exon extends GenomicEntity {
  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  public Exon(GenomicRegion region) {
    this(region, Strand.SENSE);
  }

  public Exon(GenomicRegion region, Strand strand) {
    super(GenomicType.EXON, region, strand);
  }

  public Exon(Genome genome, Chromosome chr, int start, int end) {
    super(GenomicType.EXON, genome, chr, start, end);
  }

  public Exon(Genome genome, Chromosome chr, int start, int end, Strand strand) {
    super(GenomicType.EXON, genome, chr, start, end, strand);
  }
}
