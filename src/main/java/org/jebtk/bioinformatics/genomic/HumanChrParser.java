package org.jebtk.bioinformatics.genomic;

import java.util.regex.Matcher;

public class HumanChrParser extends ChrParser {

  /**
   * Look for the numerical part of the chromosome to give it a numerical order.
   * 
   * @param name
   * @return
   */
  @Override
  public int getId(String name) {
    Matcher matcher = Chromosome.CHR_NUM_GROUP_REGEX.matcher(name);

    if (matcher.find()) {
      return Integer.parseInt(matcher.group(1));
    } else {
      char c = name.charAt(0);

      switch (c) {
      case 'X':
        return 23;
      case 'Y':
        return 24;
      default:
        // M
        return 25;
      }
    }
  }
}
