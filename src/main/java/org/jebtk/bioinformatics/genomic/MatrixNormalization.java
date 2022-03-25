package org.jebtk.bioinformatics.genomic;

import java.util.List;

import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.math.MathUtils;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.math.matrix.MatrixDimFunction;
import org.jebtk.math.matrix.utils.MatrixOperations.ColScale;

public class MatrixNormalization {

  private static class TPMRPK implements MatrixDimFunction {
    private double[] mWidthsKb;

    public TPMRPK(final double[] widths) {
      mWidthsKb = widths;
    }

    @Override
    public void apply(int index, double[] data, double[] ret) {
      MathUtils.multiply(data, mWidthsKb, ret);
    }
  }

  private static class TPMFactors implements MatrixDimFunction {
    @Override
    public void apply(int index, double[] data, double[] ret) {
      double factor = 0.0;

      for (int i = 0; i < data.length; ++i) {
        factor += data[i];
      }

      ret[index] = factor / 1000000;
    }
  }

  public static DataFrame tpm(final DataFrame m, final List<GenomicRegion> locations) {

    DataFrame ret = new DataFrame(m);

    double[] widthsKb = getWidthsKb(locations);

    int c = ret.getCols();

    ret.colApply(new TPMRPK(widthsKb));

    // Get the scaling factor for each column
    double[] factors = new double[c];
    ret.colEval(new TPMFactors(), factors);

    // Scale the values
    ret.apply(new ColScale(factors));

    return ret;
  }

  private static double[] getWidthsKb(List<GenomicRegion> locations) {
    double[] ret = new double[locations.size()];

    for (int i = 0; i < locations.size(); ++i) {
      ret[i] = locations.get(i).getLength() / 1000.0;
    }

    return ret;
  }

  public static DataFrame rpm(final DataFrame m) {
    double[] counts = m.getColumnHeader().getAnnotation("total-reads").rowToDouble(0);

    return rpm(m, counts);
  }

  public static DataFrame rpm(final DataFrame m, final List<Integer> counts) {
    return rpm(m, CollectionUtils.toDoublePrimitive(counts));
  }

  public static DataFrame rpm(final DataFrame m, final double[] counts) {
    DataFrame ret = new DataFrame(m);

    double[] factors = new double[counts.length];

    MathUtils.divide(counts, 1000000, factors);

    ret.apply(new ColScale(factors));

    return ret;
  }

  /**
   * Convert a table of counts into rpkm values.
   * 
   * @param m
   * @param counts
   * @param locations
   * @return
   */
  public static DataFrame rpkm(final DataFrame m, final List<GenomicRegion> locations) {
    DataFrame ret = rpm(m);

    double[] factors = getWidthsKb(locations);

    ret.apply(new ColScale(factors));

    return ret;
  }

  public static DataFrame rpkm(final DataFrame m, final List<Integer> counts, final List<GenomicRegion> locations) {
    DataFrame ret = rpm(m, counts);

    double[] factors = getWidthsKb(locations);

    ret.apply(new ColScale(factors));

    return ret;
  }
}
