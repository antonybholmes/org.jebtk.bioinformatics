package org.jebtk.bioinformatics.genomic;

public enum GenomicType {
  GENE, TRANSCRIPT, EXON, UTR_5P, UTR_3P, CYTOBAND, SUPER_ENHANCER, ENHANCER, PEAK, REGION;

  public static GenomicType parse(String s) {
    switch (s.toLowerCase()) {
    case "gene":
      return GENE;
    case "transcript":
      return TRANSCRIPT;
    case "CDS":
      return TRANSCRIPT;
    case "exon":
      return EXON;
    case "5p_utr":
      return UTR_5P;
    case "3p_utr":
      return UTR_3P;
    case "super_enhancer":
      return SUPER_ENHANCER;
    case "peak":
      return PEAK;
    default:
      return REGION;
    }
  }

  public static boolean le(GenomicType l1, GenomicType l2) {
    return l1.compareTo(l2) < 0;
  }

  public static boolean lt(GenomicType l1, GenomicType l2) {
    int c = l1.compareTo(l2);

    return c == 0 || c < 0;
  }
}
