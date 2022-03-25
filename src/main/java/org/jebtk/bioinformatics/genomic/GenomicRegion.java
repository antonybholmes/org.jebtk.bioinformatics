/**
 * Copyright (C) 2016, Antony Holmes
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  3. Neither the name of copyright holder nor the names of its contributors 
 *     may be used to endorse or promote products derived from this software 
 *     without specific prior written permission. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 */
package org.jebtk.bioinformatics.genomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.Formatter;
import org.jebtk.core.text.Formatter.NumberFormatter;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

import com.fasterxml.jackson.annotation.JsonGetter;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;

/**
 * Describes a region of a genome.
 *
 * @author Antony Holmes
 */
@JsonPropertyOrder({ "loc", "strand" })
public class GenomicRegion extends Region {

  private static final long serialVersionUID = 1L;

  private static final String CHR_REGEX = "(chr[a-zA-Z0-9\\-\\_]+)";

  /** The Constant GENOMIC_REGEX. */
  private static final Pattern GENOMIC_PATTERN = Pattern.compile(CHR_REGEX + ":(\\d+(?:,\\d+)*)"); // Pattern.compile("chr.+?:(?:\\d+(?:,\\d+)*)(?:-\\d+(?:,\\d+)*)?");

  public static final Pattern CHR_ONLY_PATTERN = Pattern.compile(CHR_REGEX);

  private static final Pattern GENOMIC_NUM_PATTERN = Pattern.compile("[\\-\\+\\t\\=\\: ](\\d+(?:,\\d+)*)");

  /**
   * The member chr.
   */
  @JsonIgnore
  public final Chromosome mChr;

  /**
   * The member strand.
   */
  @JsonIgnore
  public final Strand mStrand;

  // @JsonIgnore
  // public final Genome mGenome;

  // protected final String mType;

  /**
   * Instantiates a new genomic region.
   *
   * @param region the region
   */
  public GenomicRegion(GenomicRegion region) {
    this(region, region.mStrand);
  }

  public GenomicRegion(GenomicRegion region, Strand strand) {
    this(region.mChr, region.mStart, region.mEnd, strand);
  }

  public GenomicRegion(GenomicRegion region, int start, int end) {
    this(region.mChr, start, end, region.mStrand);
  }

  public GenomicRegion(Chromosome chr, int start, int end) {
    this(chr, start, end, Strand.SENSE);
  }

  public GenomicRegion(Chromosome chr, int start, int end, Strand strand) {
    super(start, end);

    mChr = chr;
    mStrand = strand;
  }

  public GenomicRegion(Genome genome, String chr, int start, int end) {
    this(ChromosomeService.getInstance().chr(genome, chr), start, end);
  }

  /**
   * Gets the strand.
   *
   * @return the strand
   */
  @JsonIgnore
  public Strand getStrand() {
    return mStrand;
  }

  /**
   * Return '+' for sense and '-' for anti-sense.
   * 
   * @return
   */
  @JsonGetter("strand")
  public char getFormattedStrand() {
    return Strand.toChar(mStrand);
  }

  /**
   * Gets the location.
   *
   * @return the location
   */
  @JsonGetter("loc")
  public String getLocation() {
    return toLocation(mChr, mStart, mEnd);
  }

  /**
   * Returns the numerical range of the location in form "<start>-<end>"
   * 
   * @return
   */
  @JsonIgnore
  public String getRange() {
    return toRange(mStart, mEnd);
  }

  /**
   * Gets the formatted location.
   *
   * @return the formatted location
   */
  @JsonIgnore
  public String getFormattedLocation() {
    return formattedLocation(mChr, mStart, mEnd);
  }

  /**
   * Returns the start considering the strand orientation. Reverse strand elements
   * return their end as their start. Relevant when comparing distances between
   * elements and closest starts.
   * 
   * @return
   */
  @JsonIgnore
  public int getStrStart() {
    if (mStrand == Strand.SENSE) {
      return mStart;
    } else {
      return mEnd;
    }
  }

  /**
   * Returns the end considering the strand orientation. Reverse strand elements
   * return their end as their start. Relevant when comparing distances between
   * elements and closest starts.
   * 
   * @return
   */
  @JsonIgnore
  public int getStrEnd() {
    if (mStrand == Strand.SENSE) {
      return mEnd;
    } else {
      return mStart;
    }
  }

  /**
   * Gets the chr.
   *
   * @return the chr
   */
  @JsonIgnore
  public Chromosome getChr() {
    return mChr;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(Region r) {
    if (r instanceof GenomicRegion) {
      GenomicRegion gr = (GenomicRegion) r;

      // int c = mChr.mGenome.compareTo(gr.mChr.mGenome); //
      // getGenome().compareTo(gr.getGenome());

      // if (c != 0) {
      // return c;
      // }

      int c = mChr.compareTo(gr.mChr);

      if (c != 0) {
        return c;
      }

      c = mStrand.compareTo(gr.mStrand);

      if (c != 0) {
        return c;
      }
    }

    // Chr are the same so look at coordinates
    return super.compareTo(r);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  @JsonIgnore
  public String toString() {
    return getLocation();
  }

  /**
   * Formatted txt.
   *
   * @return the string
   */
  public String formattedTxt() {
    StringBuilder buffer = new StringBuilder();

    try {
      formattedTxt(buffer);
      buffer.append(TextUtils.NEW_LINE);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return buffer.toString();
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.genome.Region#formattedTxt(java.lang.
   * Appendable)
   */
  @Override
  public void formattedTxt(Appendable buffer) throws IOException {
    buffer.append(mChr.toString());
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Integer.toString(mStart));
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Integer.toString(mEnd));
    buffer.append(TextUtils.NEW_LINE);
  }

  /**
   * Parses the.
   *
   * @param locations the locations
   * @return the list
   * @throws ParseException the parse exception
   */
  public static List<GenomicRegion> parse(Genome genome, List<String> locations) {
    List<GenomicRegion> ret = new ArrayList<GenomicRegion>();

    for (String location : locations) {
      GenomicRegion region = parse(genome, location);

      if (region != null) {
        ret.add(region);
      }
    }

    return ret;
  }

  /**
   * Parse a position annotation in the form (chr)n:x-(y) {.
   *
   * @param location the location
   * @return the genomic region
   */
  public static GenomicRegion parse(Genome genome, String location) {
    // System.err.println("location: " + location);

    if (Io.isEmptyLine(location)) {
      return null;
    }

    if (location.contains(TextUtils.NA)) {
      return null;
    }

    if (location.length() == 0) {
      return null;
    }

    Matcher matcher = CHR_ONLY_PATTERN.matcher(location);

    if (matcher.find()) {
      Chromosome chr = ChromosomeService.getInstance().chr(matcher.group(1));

      if (chr == null) {
        return null;
      }

      int start = 1;
      int end = 1;

      matcher = GENOMIC_NUM_PATTERN.matcher(location);

      if (matcher.find()) {
        start = TextUtils.parseInt(matcher.group(1));
      } else {
        // Just the name, return the full chromosome
        return new GenomicRegion(chr, 1, ChromosomeService.getInstance().size(genome, chr));
      }

      if (matcher.find()) {
        end = TextUtils.parseInt(matcher.group(1));
      } else {
        // If no end is specified, make the end at least the start
        end = start;
      }

      return new GenomicRegion(chr, start, end);
    } else if (isRegion(location)) {
      location = region(location);

      int start;
      int end;

      if (location.indexOf("-") != -1) {
        List<String> tokens = Splitter.on('-').text(location); // )
                                                               // .(tokens.get(1),
                                                               // '-');

        start = TextUtils.parseInt(tokens.get(0));
        end = TextUtils.parseInt(tokens.get(1));
      } else {
        // single position

        start = TextUtils.parseInt(location);
        end = start;
      }

      return GenomicRegion.create(start, end);
    } else {
      return null;
    }
  }

  public static GenomicRegion parse(Genome genome, String chr, int start, int end) {
    return new GenomicRegion(genome, chr, start, end);
  }

  /**
   * Return the overlapping region between two regions or null if they are not
   * overlapping.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the genomic region
   */
  public static GenomicRegion overlap(GenomicRegion region1, GenomicRegion region2) {
    if (!region1.mChr.equals(region2.mChr)) {
      return null;
    }

    return overlap(region1.mChr, region1.mStart, region1.mEnd, region2.mStart, region2.mEnd);
  }

  /**
   * Overlap.
   *
   * @param chr     the chr
   * @param start1  the start 1
   * @param end1    the end 1
   * @param region2 the region 2
   * @return the genomic region
   */
  public static GenomicRegion overlap(Chromosome chr, int start1, int end1, GenomicRegion region2) {
    if (!chr.equals(region2.mChr)) {
      return null;
    }

    return overlap(chr, start1, end1, region2.mStart, region2.mEnd);
  }

  /**
   * Overlap.
   *
   * @param chr    the chr
   * @param start1 the start 1
   * @param end1   the end 1
   * @param start2 the start 2
   * @param end2   the end 2
   * @return the genomic region
   */
  public static GenomicRegion overlap(Chromosome chr, int start1, int end1, int start2, int end2) {
    if (end1 < start2 || end2 < start1) {
      return null;
    }

    int start = Math.max(start1, start2);
    int end = Math.min(end1, end2);

    return new GenomicRegion(chr, start, end);
  }

  /**
   * Overlap type.
   *
   * @param region1 the region 1
   * @param region2 the region 2
   * @return the overlap type
   */
  public static OverlapType overlapType(GenomicRegion region1, GenomicRegion region2) {
    if (!region1.mChr.equals(region2.mChr)) {
      return OverlapType.NONE;
    }

    return overlapType(region1.mChr, region1.mStart, region1.mEnd, region2.mStart, region2.mEnd);
  }

  /**
   * Overlap type.
   *
   * @param chr    the chr
   * @param start1 the start 1
   * @param end1   the end 1
   * @param start2 the start 2
   * @param end2   the end 2
   * @return the overlap type
   */
  public static OverlapType overlapType(Chromosome chr, int start1, int end1, int start2, int end2) {
    if (end1 < start2 || end2 < start1) {
      return OverlapType.NONE;
    }

    if (start1 >= start2 && end1 <= end2) {
      return OverlapType.WITHIN;
    } else if (start1 < start2 && end1 > end2) {
      return OverlapType.COMPLETE;
    } else {
      return OverlapType.PARTIAL;
    }
  }

  /**
   * Overlaps.
   *
   * @param chr    the chr
   * @param start1 the start 1
   * @param end1   the end 1
   * @param start2 the start 2
   * @param end2   the end 2
   * @return true, if successful
   */
  public static boolean overlaps(Chromosome chr, int start1, int end1, int start2, int end2) {
    return (start1 >= start2 && start1 < end2) || (end1 >= start2 && end1 < end2) || (start2 >= start1 && start2 < end1)
        || (end2 >= start1 && end2 < end1);
  }

  /**
   * Use only if start1 is before start2.
   *
   * @param start1 the start 1
   * @param end1   the end 1
   * @param start2 the start 2
   * @param end2   the end 2
   * @return true, if successful
   */
  public static boolean overlapOrdereds(int start1, int end1, int start2, int end2) {
    return end1 >= start2;
  }

  /**
   * Returns true if region2 overlaps with region1.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return true, if successful
   */
  public static boolean overlaps(GenomicRegion region1, GenomicRegion region2) {
    return overlap(region1, region2) != null;
  }

  /**
   * Returns the mid point of a region.
   *
   * @param region the region
   * @return the int
   */
  public static int mid(GenomicRegion region) {
    return mid(region.mStart, region.mEnd); // new GenomicRegion(region.mChr,
                                            // mid, mid);
  }

  /**
   * Mid.
   *
   * @param start the start
   * @param end   the end
   * @return the int
   */
  public static int mid(int start, int end) {
    return (start + end) / 2; // new GenomicRegion(region.mChr, mid, mid);
  }

  /**
   * Return the mid point of a region i.e. the start and end will be the same mid
   * point of the genomic region.
   *
   * @param region the region
   * @return the genomic region
   */
  public static GenomicRegion midRegion(GenomicRegion region) {
    if (region == null) {
      return null;
    }

    int mid = mid(region);

    return new GenomicRegion(region.mChr, mid, mid);
  }

  /**
   * Return the distance between the midpoint of region1 and region2. If region1
   * is greater than region2, a positive value will be returned.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the int
   */
  public static int midDist(GenomicRegion region1, GenomicRegion region2) {
    return mid(region1) - mid(region2);
  }

  /**
   * Mid abs dist.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the int
   */
  public static int midAbsDist(GenomicRegion region1, GenomicRegion region2) {
    return Math.abs(midDist(region1, region2));
  }

  /**
   * Returns true if the string is of the form chrX:A-B where X is either a number
   * or letter and A and B are numbers.
   *
   * @param text the text
   * @return true, if is region
   */
  public static boolean isGenomicRegion(String text) {
    if (text == null) {
      return false;
    }

    return GENOMIC_PATTERN.matcher(text).find(); // value.matches("chr.+?:\\d+(-\\d+)?.*");
  }

  public static String genomicRegion(String text) {
    if (text == null) {
      return TextUtils.EMPTY_STRING;
    }

    Matcher matcher = GENOMIC_PATTERN.matcher(text);

    if (matcher.find()) {
      return matcher.group(1);
    } else {
      return TextUtils.EMPTY_STRING;
    }
  }

  /**
   * Returns true if the value appears to be a chromosome.
   *
   * @param value the value
   * @return true, if is chr
   */
  public static boolean isChr(String value) {
    if (value == null) {
      return false;
    }

    return value.startsWith("chr");
  }

  /**
   * Returns true if a point is within a region.
   *
   * @param start  the start
   * @param region the region
   * @return true, if successful
   */
  public static boolean within(int start, GenomicRegion region) {
    return start >= region.mStart && start <= region.mEnd;
  }

  /**
   * Extend.
   *
   * @param regions     the regions
   * @param startOffset the start offset
   * @param endOffset   the end offset
   * @return the list
   */
  public static List<GenomicRegion> extend(List<GenomicRegion> regions, int startOffset, int endOffset) {
    List<GenomicRegion> ret = new ArrayList<GenomicRegion>();

    for (GenomicRegion region : regions) {
      ret.add(extend(region, startOffset, endOffset));
    }

    return ret;
  }

  /**
   * Extend.
   *
   * @param region      the region
   * @param startOffset the start offset
   * @param endOffset   the end offset
   * @return the genomic region
   */
  public static GenomicRegion extend(GenomicRegion region, int startOffset, int endOffset) {
    return new GenomicRegion(region.mChr, region.mStart - Math.abs(startOffset), region.mEnd + Math.abs(endOffset));
  }

  /**
   * Extend around a position.
   *
   * @param chr         the chr
   * @param start       the start
   * @param startOffset the start offset
   * @param endOffset   the end offset
   * @return the genomic region
   */
  public static GenomicRegion extend(Genome genome, Chromosome chr, int start, int startOffset, int endOffset) {
    return new GenomicRegion(chr, start - Math.abs(startOffset), start + Math.abs(endOffset));
  }

  /**
   * Return the sequences for a set of regions.
   *
   * @param genome         the genome
   * @param regions        the regions
   * @param genomeAssembly the genome assembly
   * @return the sequences
   * @throws IOException    Signals that an I/O exception has occurred.
   * @throws ParseException the parse exception
   */
  public static List<SequenceRegion> getSequences(Genome genome, List<GenomicRegion> regions,
      SequenceReader genomeAssembly) throws IOException, ParseException {
    return getSequences(genome, regions, true, RepeatMaskType.UPPERCASE, genomeAssembly);
  }

  /**
   * Gets the sequences.
   *
   * @param genome         the genome
   * @param regions        the regions
   * @param displayUpper   the display upper
   * @param repeatMaskType the repeat mask type
   * @param genomeAssembly the genome assembly
   * @return the sequences
   * @throws IOException    Signals that an I/O exception has occurred.
   * @throws ParseException the parse exception
   */
  public static List<SequenceRegion> getSequences(Genome genome, Collection<GenomicRegion> regions,
      boolean displayUpper, RepeatMaskType repeatMaskType, SequenceReader genomeAssembly)
      throws IOException, ParseException {
    List<SequenceRegion> sequences = new ArrayList<SequenceRegion>();

    for (GenomicRegion region : regions) {
      System.err.println("Getting sequence for region " + region.getLocation() + "...");

      sequences.add(genomeAssembly.getSequence(genome, region, displayUpper, repeatMaskType));
    }

    return sequences;
  }

  /**
   * Sort regions by start and group by chromosome.
   *
   * @param <T>     the generic type
   * @param regions the regions
   * @return the map
   */
  public static <T extends GenomicRegion> Map<Chromosome, List<T>> sortByStart(Iterable<T> regions) {
    Map<Chromosome, Map<Integer, T>> map = new HashMap<Chromosome, Map<Integer, T>>();

    for (T region : regions) {
      if (!map.containsKey(region.mChr)) {
        map.put(region.mChr, new TreeMap<Integer, T>());
      }

      map.get(region.mChr).put(region.mStart, region);
    }

    Map<Chromosome, List<T>> ret = new HashMap<Chromosome, List<T>>();

    for (Chromosome chr : map.keySet()) {
      ret.put(chr, new ArrayList<T>());

      for (int start : map.get(chr).keySet()) {
        ret.get(chr).add(map.get(chr).get(start));
      }
    }

    return ret;
  }

  /**
   * Sorts the regions on the assumption they are on the same chromosome. Thus it
   * returns a list with the regions sorted by start.
   *
   * @param <T>     the generic type
   * @param regions the regions
   * @return the list
   */
  public static <T extends GenomicRegion> List<T> sortByStartSingleChr(List<T> regions) {
    Map<Integer, T> map = new TreeMap<Integer, T>();

    for (T region : regions) {
      map.put(region.mStart, region);
    }

    List<T> ret = new ArrayList<T>();

    for (int start : map.keySet()) {
      ret.add(map.get(start));
    }

    return ret;
  }

  /**
   * Sort single chr.
   *
   * @param <T>     the generic type
   * @param regions the regions
   * @return the list
   */
  public static <T extends GenomicRegion> List<T> sortSingleChr(List<T> regions) {
    Map<Integer, T> map = new TreeMap<Integer, T>();

    for (T region : regions) {
      map.put(region.mStart, region);
      map.put(region.mEnd, region);
    }

    List<T> ret = new ArrayList<T>();

    for (int start : map.keySet()) {
      ret.add(map.get(start));
    }

    return ret;
  }

  /**
   * Shift a region by a number of bases.
   *
   * @param region the region
   * @param shift  the shift
   * @param sizes  the sizes
   * @return the genomic region
   */
  public static GenomicRegion shift(Genome genome, GenomicRegion region, int shift) {

    return shift(genome, region, shift, 0);

  }

  public static GenomicRegion shift(Genome genome, GenomicRegion region, int shift, int minSep) {

    int size = ChromosomeService.getInstance().size(genome, region.mChr);

    // bound the positions so they dont exceed the chromosome bounds
    int start = Math.min(size, Math.max(1, region.mStart + shift));
    int end = Math.min(size, Math.max(start + minSep, region.mEnd + shift));

    return new GenomicRegion(region.mChr, start, end);

  }

  public static GenomicRegion add(GenomicRegion region, int shift) {

    return add(region, shift, 0);
  }

  public static GenomicRegion add(GenomicRegion region, int shift, int minSep) {

    // bound the positions so they dont exceed the chromosome bounds
    int start = Math.max(1, region.mStart + shift);
    int end = Math.max(start + minSep, region.mEnd + shift);

    return new GenomicRegion(region.mChr, start, end);

  }

  /**
   * Returns true if region1 is within region2 (does not extend outside boundary).
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return true, if successful
   */
  public static boolean within(GenomicRegion region1, GenomicRegion region2) {
    if (region1 == null || region2 == null) {
      return false;
    }

    return region1.mChr.equals(region2.mChr) && region1.mStart >= region2.mStart && region1.mEnd <= region2.mEnd;
  }

  /**
   * Returns the closest absolute distance between two regions.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the int
   */
  public static int minAbsDist(GenomicRegion region1, GenomicRegion region2) {
    return Math.min(Math.abs(minDist(region1, region2)), Math.abs(minDist(region2, region1)));
  }

  /**
   * Min dist.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the int
   */
  public static int minDist(GenomicRegion region1, GenomicRegion region2) {
    return region2.mStart - region1.mEnd;
  }

  /**
   * Returns the width of region1 as a percentage of region2.
   *
   * @param region1 the region1
   * @param region2 the region2
   * @return the double
   */
  public static double p(GenomicRegion region1, GenomicRegion region2) {
    return (double) region1.mLength / (double) region2.mLength;
  }

  /**
   * To location.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the string
   */
  public static String toLocation(Chromosome chr, int start, int end) {
    return chr + ":" + toRange(start, end);
  }

  /**
   * Formatted location.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the string
   */
  public static String formattedLocation(Chromosome chr, int start, int end) {
    NumberFormatter f = Formatter.number();

    return chr + ":" + f.format(start) + "-" + f.format(end);
  }

  /**
   * Creates a special kind of genomic region with no chromosome. Should be used
   * in situations where the chromosome is irrelevent, but the method call expects
   * it.
   *
   * @param start the start
   * @param end   the end
   * @return A genomic region with an invalid chromosome.
   */
  public static GenomicRegion create(int start, int end) {
    return create(Chromosome.NO_CHR, start, end);
  }

  /*
   * public static GenomicRegion create(String genome, String region) { int s = 0;
   * int e = region.indexOf(':');
   * 
   * String chr = region.substring(0, e);
   * 
   * s = e + 1;
   * 
   * e = region.lastIndexOf('-');
   * 
   * String start; String end;
   * 
   * if (e != -1) { start = region.substring(s, e);
   * 
   * s = e + 1;
   * 
   * end = region.substring(s);
   * 
   * } else { start = region.substring(s); end = start; }
   * 
   * 
   * return create(genome, chr, start, end); }
   */

  public static GenomicRegion create(String chr, String start, String end) {
    return create(chr, TextUtils.parseInt(start), TextUtils.parseInt(end));
  }

  public static GenomicRegion create(String chr, int start, int end) {
    return create(ChromosomeService.getInstance().chr(chr), start, end);
  }

  /**
   * Creates the.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the genomic region
   */
  public static GenomicRegion create(Chromosome chr, int start, int end) {
    return create(chr, start, end, Strand.SENSE);
  }

  /**
   * Creates the.
   *
   * @param chr    the chr
   * @param start  the start
   * @param end    the end
   * @param strand the strand
   * @return the genomic region
   */
  public static GenomicRegion create(Chromosome chr, int start, int end, Strand strand) {
    return new GenomicRegion(chr, start, end, strand);
  }

  /**
   * Center regions.
   *
   * @param regions the regions
   * @return the list
   */
  public static List<GenomicRegion> center(List<GenomicRegion> regions) {
    List<GenomicRegion> ret = new ArrayList<GenomicRegion>(regions.size());

    for (GenomicRegion region : regions) {
      ret.add(center(region));
    }

    return ret;
  }

  /**
   * Center.
   *
   * @param region the region
   * @return the genomic region
   */
  public static GenomicRegion center(GenomicRegion region) {
    int mid = mid(region);

    return new GenomicRegion(region.mChr, mid, mid);
  }

  public static GenomicRegion randomRegion(Genome genome, int length) throws IOException {
    Random rand = new Random();

    Chromosome chr = ChromosomeService.getInstance().randChr(genome);

    int start = rand.nextInt(ChromosomeService.getInstance().size(genome, chr) - length);
    int end = start + length - 1;

    return new GenomicRegion(chr, start, end);
  }

  /**
   * Returns a copy of the region with the strand changed to its oppoite.
   * 
   * @param region
   * @return
   */
  public static GenomicRegion oppositeStrand(GenomicRegion region) {
    return new GenomicRegion(region.mChr, region.mStart, region.mEnd, Strand.oppositeStrand(region.mStrand));
  }

  public static GenomicRegion reverseCompliment(GenomicRegion region) {
    return new GenomicRegion(region.mChr, region.mEnd, region.mStart, Strand.oppositeStrand(region.mStrand));
  }

  public static List<GenomicRegion> parseRegions(Genome genome, Path file) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    List<GenomicRegion> ret = new ArrayList<GenomicRegion>();

    String line;

    try {
      reader.readLine();

      while ((line = reader.readLine()) != null) {
        ret.add(parse(genome, line));
      }
    } finally {
      reader.close();
    }

    return ret;
  }
}
