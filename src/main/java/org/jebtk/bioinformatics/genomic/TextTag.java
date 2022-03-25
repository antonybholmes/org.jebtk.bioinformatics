package org.jebtk.bioinformatics.genomic;

public class TextTag extends Tag {
  private String mValue;

  public TextTag(String value) {
    mValue = value;
  }

  @Override
  public String toString() {
    return mValue;
  }

  @Override
  public int toInt() {
    return Integer.parseInt(mValue);
  }

  @Override
  public double toDouble() {
    return Double.parseDouble(mValue);
  }

  @Override
  public TagType getType() {
    return TagType.TEXT;
  }

  @Override
  protected int compareValue(Tag p) {
    return mValue.compareTo(p.toString());
  }
}
