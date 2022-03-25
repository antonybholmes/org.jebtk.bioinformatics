package org.jebtk.bioinformatics.annotation;

/**
 * The Class Species.
 */
public class Genome extends Type {

  /**
   * Instantiates a new species.
   *
   * @param name the name
   */
  public Genome(String name) {
    this(-1, name);
  }

  /**
   * Instantiates a new species.
   *
   * @param id   the id
   * @param name the name
   */
  public Genome(int id, String name) {
    super(id, name);
  }
}
