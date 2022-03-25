package org.jebtk.bioinformatics.genomic;

public class DoubleTag extends Tag {
  private double mValue;

  public DoubleTag(double value) {
    mValue = value;
  }

  @Override
  public int toInt() {
    return (int) mValue;
  }

  @Override
  public double toDouble() {
    return mValue;
  }

  @Override
  public String toString() {
    return Double.toString(mValue);
  }

  @Override
  public TagType getType() {
    return TagType.DOUBLE;
  }

  @Override
  protected int compareValue(Tag p) {
    double v = p.toDouble();

    if (mValue > v) {
      return 1;
    } else if (mValue < v) {
      return -1;
    } else {
      return 0;
    }
  }
}
