package org.jebtk.bioinformatics.genomic;

public class IntTag extends Tag {
  private int mValue;

  public IntTag(int value) {
    mValue = value;
  }

  @Override
  public int toInt() {
    return mValue;
  }

  @Override
  public double toDouble() {
    return mValue;
  }

  @Override
  public String toString() {
    return Integer.toString(mValue);
  }

  @Override
  public TagType getType() {
    return TagType.INT;
  }

  @Override
  protected int compareValue(Tag p) {
    int v = p.toInt();

    if (mValue > v) {
      return 1;
    } else if (mValue < v) {
      return -1;
    } else {
      return 0;
    }
  }
}
