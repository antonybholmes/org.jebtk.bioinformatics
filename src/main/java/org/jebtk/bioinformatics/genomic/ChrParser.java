package org.jebtk.bioinformatics.genomic;

import java.util.regex.Matcher;

public class ChrParser {

  /**
   * Look for the numerical part of the chromosome to give it a numerical order.
   * 
   * Numerical chrosomes number 0 to 2^16 (exclusive). Lettered chrosomes are
   * numbered 2^16 to 2^24. Chromosomes containng '_'
   * 
   * @param name
   * @return
   */
  public int getId(String name) {
    Matcher matcher = Chromosome.CHR_NUM_GROUP_REGEX.matcher(name);

    int ret;

    if (matcher.find()) {
      ret = Short.parseShort(matcher.group(1));
    } else {
      // Assume a letter so subtract 65 (so A becomes 1, B 2 etc) and then
      // Shift 16 bits,
      ret = (name.charAt(0) - 64) << 16;
    }

    return ret;
  }
}
