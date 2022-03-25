package org.jebtk.bioinformatics.genomic;

public enum TagType {
  TEXT, INT, DOUBLE;

  public static byte byteRep(TagType type) {
    switch (type) {
    case INT:
      return 1;
    case DOUBLE:
      return 2;
    default:
      // String
      return 0;
    }
  }

  public static TagType parse(int b) {
    switch (b) {
    case 2:
      return DOUBLE;
    case 1:
      return INT;
    default:
      return TEXT;
    }
  }
}
