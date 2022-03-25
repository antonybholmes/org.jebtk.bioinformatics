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
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.bioinformatics.genomic.GenomeService;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Describes a BED file.
 * 
 * @author Antony Holmes
 *
 */
public class BedGraph extends UCSCTrack {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant HEIGHT_PATTERN.
   */
  public static final Pattern HEIGHT_PATTERN = Pattern.compile("maxHeightPixels=(\\d+):(\\d+):(\\d+)");

  /**
   * The constant LOG.
   */
  private static final Logger LOG = LoggerFactory.getLogger(BedGraph.class);

  /**
   * The constant DEFAULT_BEDGRAPH_COLOR.
   */
  protected static final Color DEFAULT_BEDGRAPH_COLOR = Color.RED;

  /**
   * The constant TRACK_TYPE.
   */
  private static final String TRACK_TYPE = "bedGraph";

  public BedGraph(String name) {
    this(name, name);
  }

  public BedGraph(String name, String description) {
    this(name, description, DEFAULT_BEDGRAPH_COLOR, DEFAULT_HEIGHT);
  }

  /**
   * Instantiates a new bed graph.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   */
  public BedGraph(String name, String description, Color color) {
    this(name, description, color, DEFAULT_HEIGHT);
  }

  /**
   * Instantiates a new bed graph.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @param height      the height
   */
  public BedGraph(String name, String description, Color color, int height) {
    super(name, description, color, TRACK_TYPE);

    mHeight = height;
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.external.ucsc.UCSCTrack#bufferHeader(
   * int, java.lang.Appendable)
   */
  @Override
  public void bufferHeader(int priority, Appendable buffer) throws IOException {
    super.bufferHeader(priority, buffer);

    buffer.append(" maxHeightPixels=").append("128:").append(Integer.toString(mHeight)).append(":1");
    buffer.append(" visibility=full autoScale=on alwaysZero=on");
  }

  /**
   * Parses the.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<BedGraph> parse(Path file) throws IOException {
    LOG.info("Parsing BedGraph file {}...", file);

    BufferedReader reader = FileUtils.newBufferedReader(file);

    BedGraph bedgraph = new BedGraph(PathUtils.getNameNoExt(file));

    String line;
    Matcher matcher;

    List<BedGraph> tracks = new ArrayList<BedGraph>();

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

          Color color = DEFAULT_BEDGRAPH_COLOR;

          if (matcher.find()) {
            color = parseColor(matcher);
          }

          // Check if we can find the height
          int height = DEFAULT_HEIGHT;

          matcher = HEIGHT_PATTERN.matcher(line);

          if (matcher.find()) {
            height = Integer.parseInt(matcher.group(2));
          }

          bedgraph = new BedGraph(name, description, color, height);

          tracks.add(bedgraph);
        } else {
          BedGraphElement region = BedGraphElement.parse(GenomeService.getInstance().guessGenome(file), line);

          bedgraph.add(region);
        }
      }
    } finally {
      reader.close();
    }

    if (tracks.size() == 0) {
      tracks.add(bedgraph);
    }

    LOG.info("BED {} ({} peaks).", bedgraph.getName(), bedgraph.size());

    return tracks;
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
      // Skip header
      reader.readLine();

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        if (isTrackLine(line)) {
          break;
        }

        ++r;
      }
    } finally {
      reader.close();
    }

    DataFrame ret = DataFrame.createDataFrame(r, 4);

    ret.setColumnName(0, "Chr");
    ret.setColumnName(1, "Start");
    ret.setColumnName(2, "End");
    ret.setColumnName(3, "Value");

    reader = FileUtils.newBufferedReader(file);

    r = 0;

    try {
      // Skip header
      reader.readLine();

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        if (isTrackLine(line)) {
          break;
        }

        List<String> tokens = TextUtils.tabSplit(line);

        ret.set(r, 0, tokens.get(0));
        ret.set(r, 1, Integer.parseInt(tokens.get(1)));
        ret.set(r, 2, Integer.parseInt(tokens.get(2)));
        ret.set(r, 3, Double.parseDouble(tokens.get(3)));

        ++r;

      }
    } finally {
      reader.close();
    }

    return ret;
  }

  public BedGraph getBedGraph(GenomicRegion region) {
    BedGraph ret = new BedGraph(mName, mDescription, mColor, mHeight);

    ret.addAll(find(region));

    return ret;
  }

}
