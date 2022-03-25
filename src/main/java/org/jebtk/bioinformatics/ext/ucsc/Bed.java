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
package org.jebtk.bioinformatics.ext.ucsc;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomeService;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.math.matrix.MixedMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Describes a BED file (a list of genomic positions with a track).
 * 
 * @author Antony Holmes
 *
 */
public class Bed extends UCSCTrack {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(Bed.class);

  /**
   * The constant DEFAULT_BED_COLOR.
   */
  public static final Color DEFAULT_BED_COLOR = Color.BLUE;

  /**
   * The constant BED_TRACK_TYPE.
   */
  public static final String BED_TRACK_TYPE = "bed";

  /**
   * Instantiates a new bed.
   *
   * @param name the name
   */
  public Bed(String name) {
    this(name, name);
  }

  /**
   * Instantiates a new bed.
   *
   * @param name        the name
   * @param description the description
   */
  public Bed(String name, String description) {
    this(name, description, DEFAULT_BED_COLOR);
  }

  /**
   * Instantiates a new bed.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   */
  public Bed(String name, String description, Color color) {
    super(name, description, color, BED_TRACK_TYPE);
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.external.ucsc.UCSCTrack#setColor(java.
   * awt .Color)
   */
  @Override
  public void setColor(Color color) {
    // Override the color preferences for all segments
    for (Entry<Chromosome, List<GenomicElement>> item : mRegions) {
      for (GenomicElement r : item.getValue()) {
        ((BedElement) r).setColor(color);
      }
    }

    super.setColor(color);
  }

  /**
   * Parses the track.
   *
   * @param file the file
   * @return the UCSC track
   * @throws IOException    Signals that an I/O exception has occurred.
   * @throws ParseException the parse exception
   */
  public static UCSCTrack parseTrack(GenomicType type, Path file) throws IOException {
    return parseTracks(type, file).get(0);
  }

  /**
   * Parses the.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<UCSCTrack> parseTracks(GenomicType type, Path file) throws IOException {
    LOG.info("Parsing BED file {}...", file);

    BufferedReader reader = FileUtils.newBufferedReader(file);

    return parseTracks(type, GenomeService.getInstance().guessGenome(file), getName(file), reader);
  }

  /**
   * Parses the.
   *
   * @param defaultName the default name
   * @param reader      the reader
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<UCSCTrack> parseTracks(GenomicType type, Genome genome, String defaultName, BufferedReader reader)
      throws IOException {
    Bed bed = null;

    String line;
    Matcher matcher;

    List<UCSCTrack> tracks = new ArrayList<UCSCTrack>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        matcher = NAME_PATTERN.matcher(line);

        if (matcher.find()) {

          String name = matcher.group(1);

          String description = name;

          matcher = DESCRIPTION_PATTERN.matcher(line);

          if (matcher.find()) {
            description = matcher.group(1);
          }

          matcher = COLOR_PATTERN.matcher(line);

          Color color = DEFAULT_BED_COLOR;

          if (matcher.find()) {
            color = parseColor(matcher);
          }

          bed = new Bed(name, description, color);

          tracks.add(bed);
        } else {
          if (bed == null) {
            bed = new Bed(defaultName);
            tracks.add(bed);
          }

          GenomicElement region = BedElement.parse(type, genome, line);

          if (region != null) {
            bed.add(region);
          }
        }
      }
    } finally {
      reader.close();
    }

    LOG.info("BED {} ({} peaks).", bed.getName(), bed.size());

    return tracks;
  }

  /**
   * Creates the bed from track line.
   *
   * @param line the line
   * @return the bed
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Bed createBedFromTrackLine(String line) throws IOException {
    String name = null;
    String description = name;
    Color color = DEFAULT_BED_COLOR;

    // see if the file has a header

    Matcher matcher;

    matcher = NAME_PATTERN.matcher(line);

    matcher.find();

    name = matcher.group(1);

    matcher = DESCRIPTION_PATTERN.matcher(line);

    if (matcher.find()) {
      description = matcher.group(1);
    }

    matcher = COLOR_PATTERN.matcher(line);

    if (matcher.find()) {
      color = parseColor(matcher);
    }

    matcher = NAME_PATTERN.matcher(line);

    if (matcher.find()) {
      name = matcher.group(1);
    }

    matcher = DESCRIPTION_PATTERN.matcher(line);

    if (matcher.find()) {
      description = matcher.group(1);
    }

    Bed bed = new Bed(name, description, color);

    return bed;
  }

  /**
   * Parses the bed graph.
   *
   * @param file the file
   * @return the bed
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Bed parseBedGraph(GenomicType type, Path file) throws IOException {
    return parseBedGraphs(type, file).get(0);
  }

  /**
   * Parses the bed graphs.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Bed> parseBedGraphs(GenomicType type, Path file) throws IOException {
    LOG.info("Parsing Bedgraph as BED file {}...", file);

    BufferedReader reader = Files.newBufferedReader(file, StandardCharsets.UTF_8);

    Bed bed = null;

    String line;

    List<Bed> beds = new ArrayList<Bed>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        if (isTrackLine(line)) {
          bed = createBedFromTrackLine(line);

          beds.add(bed);
        } else {
          GenomicElement region = BedElement.parse(type, GenomeService.getInstance().guessGenome(file), line);
          bed.add(region);
        }
      }
    } finally {
      reader.close();
    }

    return beds;
  }

  /**
   * Write.
   *
   * @param regions the regions
   * @param file    the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void writeBed(List<GenomicRegion> regions, Path file) throws IOException {
    writeBed(PathUtils.getName(file), regions, file);
  }

  /**
   * Write.
   *
   * @param name    the name
   * @param regions the regions
   * @param file    the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void writeBed(String name, List<GenomicRegion> regions, Path file) throws IOException {
    writeBed(name, name, regions, file);
  }

  /**
   * Write.
   *
   * @param name        the name
   * @param description the description
   * @param regions     the regions
   * @param file        the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void writeBed(String name, String description, List<GenomicRegion> regions, Path file)
      throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    try {
      writer.write(getHeader(BED_TRACK_TYPE, name, description, DEFAULT_BED_COLOR));

      writer.newLine();

      for (GenomicRegion region : regions) {
        writer.write(region.getChr().toString());
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(Integer.toString(region.getStart() - 1));
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(Integer.toString(region.getEnd()));
        writer.write(TextUtils.TAB_DELIMITER);
        writer.write(region.getLocation());
        writer.newLine();
      }
    } finally {
      writer.close();
    }
  }

  /**
   * Gets the header.
   *
   * @param name        the name
   * @param description the description
   * @return the header
   */
  public static String getHeader(String name, String description) {
    return getHeader(name, description, DEFAULT_BED_COLOR);
  }

  /**
   * Gets the header.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @return the header
   */
  public static String getHeader(String name, String description, Color color) {
    return getHeader(name, description, color, 1);
  }

  /**
   * Gets the header.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @param priority    the priority
   * @return the header
   */
  public static String getHeader(String name, String description, Color color, int priority) {
    return getHeader(BED_TRACK_TYPE, name, description, color, priority);
  }

  /**
   * Create a BED line from a region.
   *
   * @param region the region
   * @return the string
   */
  public static String toString(GenomicRegion region) {
    return toString(region, region.getLocation());
  }

  /**
   * Create a BED line from a region with a name.
   *
   * @param region the region
   * @param name   the name
   * @return the string
   */
  public static String toString(GenomicRegion region, String name) {
    StringBuilder buffer = new StringBuilder();

    buffer.append(region.getChr().toString());
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Integer.toString(region.getStart()));
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(Integer.toString(region.getEnd()));
    buffer.append(TextUtils.TAB_DELIMITER);
    buffer.append(name);

    return buffer.toString();
  }

  /**
   * To matrix.
   *
   * @param file the file
   * @return the annotation matrix
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static DataFrame toMatrix(Path file) throws IOException {
    String line;

    BufferedReader reader = FileUtils.newBufferedReader(file);

    int r = 0;

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        if (isTrackLine(line)) {
          continue;
        }

        ++r;
      }
    } finally {
      reader.close();
    }

    DataFrame ret = new DataFrame(new MixedMatrix(r, 4));

    ret.setColumnName(0, "Chr");
    ret.setColumnName(1, "Start");
    ret.setColumnName(2, "End");
    ret.setColumnName(3, "Name");

    reader = FileUtils.newBufferedReader(file);

    r = 0;

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        if (isTrackLine(line)) {
          continue;
        }

        List<String> tokens = TextUtils.tabSplit(line);

        ret.set(r, 0, tokens.get(0));
        // UCSC convention
        ret.set(r, 1, Integer.parseInt(tokens.get(1)) + 1);
        ret.set(r, 2, Integer.parseInt(tokens.get(2)));

        if (tokens.size() > 3) {
          ret.set(r, 3, tokens.get(3));
        }

        ++r;

      }
    } finally {
      reader.close();
    }

    return ret;
  }

  /**
   * Creates the.
   *
   * @param name    the name
   * @param regions the regions
   * @return the bed
   */
  public static Bed create(GenomicType type, String name, Collection<GenomicRegion> regions) {

    Bed bed = new Bed(name);

    for (GenomicRegion region : regions) {
      bed.add(BedElement.create(type, region));
    }

    return bed;
  }
}
