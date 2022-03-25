package org.jebtk.bioinformatics.genomic;

public interface ChromosomeReader extends Iterable<Chromosome> {
  public Genome getGenome();

  int size(Chromosome chr);

  public Chromosome chr(String chr);

  public Chromosome randChr();
}
