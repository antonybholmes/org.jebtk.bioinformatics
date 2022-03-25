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
import java.nio.file.Path;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicElementsMap;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.core.NameGetter;
import org.jebtk.core.Resources;
import org.jebtk.core.SizeGetter;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.event.ChangeListeners;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.TextUtils;

/**
 * The class UCSCTrack.
 */
public class UCSCTrack extends ChangeListeners
    implements NameGetter, SizeGetter, Iterable<Entry<Chromosome, List<GenomicElement>>> {

  /**
   * The constant serialVersionUID.
   */
  private static final long serialVersionUID = 1L;

  /**
   * The constant ATTRIBUTE_PATTERN.
   */
  public static final Pattern ATTRIBUTE_PATTERN = Pattern.compile("([^ ]+)=([\"'].+?[\"']|[^ ]+)");

  /**
   * The constant NAME_PATTERN.
   */
  public static final Pattern NAME_PATTERN = Pattern.compile("name=[\"'](.+?)[\"']");

  /**
   * The constant DESCRIPTION_PATTERN.
   */
  public static final Pattern DESCRIPTION_PATTERN = Pattern.compile("description=[\"'](.+?)[\"']");

  /**
   * The constant COLOR_PATTERN.
   */
  public static final Pattern COLOR_PATTERN = Pattern.compile("(\\d+),(\\d+),(\\d+)");

  /**
   * The constant TRACK_PREFIX.
   */
  protected static final String TRACK_PREFIX = "track";

  /**
   * The member name.
   */
  protected String mName;

  /**
   * The member description.
   */
  protected String mDescription;

  /**
   * The member color.
   */
  protected Color mColor;

  /**
   * The member type.
   */
  private String mType;

  /**
   * The member chr regions.
   */
  protected GenomicElementsMap mRegions = new GenomicElementsMap();

  /**
   * The constant DEFAULT_HEIGHT.
   */
  public static final int DEFAULT_HEIGHT = 128;

  /**
   * The member height.
   */
  protected int mHeight = DEFAULT_HEIGHT;

  /**
   * Instantiates a new UCSC track.
   *
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @param type        the type
   */
  public UCSCTrack(String name, String description, Color color, String type) {
    mName = name;
    mDescription = description;
    mColor = color;
    mType = type;
  }

  /**
   * Sets the name.
   *
   * @param name the new name
   */
  public void setName(String name) {
    mName = name;

    fireChanged();
  }

  /**
   * Sets the description.
   *
   * @param description the new description
   */
  public void setDescription(String description) {
    mDescription = description;

    fireChanged();
  }

  /**
   * Sets the height.
   *
   * @param height the new height
   */
  public void setHeight(int height) {
    mHeight = height;

    fireChanged();
  }

  /**
   * Gets the height.
   *
   * @return the height
   */
  public int getHeight() {
    return mHeight;
  }

  /**
   * Sets the color.
   *
   * @param color the new color
   */
  public void setColor(Color color) {
    mColor = color;

    fireChanged();
  }

  /**
   * Gets the type.
   *
   * @return the type
   */
  public String getType() {
    return mType;
  }

  /**
   * The name of the bed file.
   *
   * @return the name
   */
  public String getName() {
    return mName;
  }

  /**
   * The bed description.
   *
   * @return the description
   */
  public String getDescription() {
    return mDescription;
  }

  /**
   * Returns the track color.
   *
   * @return the color
   */
  public Color getColor() {
    return mColor;
  }

  public UCSCTrack add(GenomicElement e) {
    mRegions.add(e);

    return this;
  }

  public UCSCTrack addAll(Collection<GenomicElement> elements) {
    mRegions.addAll(elements);

    return this;
  }

  /**
   * Return the regions associated with this track.
   *
   * @return the regions
   */
  public List<GenomicElement> getElements() {
    return mRegions.getElements();
  }

  public List<GenomicElement> getElements(Chromosome chr) {
    return mRegions.getElements(chr);
  }

  public List<GenomicElement> find(GenomicRegion region) {
    return mRegions.find(region);
  }

  @Override
  public Iterator<Entry<Chromosome, List<GenomicElement>>> iterator() {
    return mRegions.iterator();
  }

  @Override
  public int size() {
    return mRegions.size();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return toString(this, 1);
  }

  /**
   * Convert track to another type.
   *
   * @param track    the track
   * @param priority the priority
   * @return the string
   */
  public String toString(UCSCTrack track, int priority) {

    StringBuilder buffer = new StringBuilder();

    try {
      toBuffer(track, priority, buffer);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return buffer.toString();
  }

  /**
   * To buffer.
   *
   * @param track    the track
   * @param priority the priority
   * @param buffer   the buffer
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public void toBuffer(UCSCTrack track, int priority, Appendable buffer) throws IOException {
    bufferHeader(priority, buffer);

    buffer.append(TextUtils.NEW_LINE);

    for (Entry<Chromosome, List<GenomicElement>> item : mRegions) {
      for (GenomicElement region : item.getValue()) {
        region.formattedTxt(buffer);
        buffer.append(TextUtils.NEW_LINE);
      }
    }
  }

  /**
   * Gets the header.
   *
   * @param priority the priority
   * @return the header
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public String getHeader(int priority) throws IOException {
    StringBuilder buffer = new StringBuilder();

    bufferHeader(priority, buffer);

    return buffer.toString();
  }

  /**
   * Buffer header.
   *
   * @param priority the priority
   * @param buffer   the buffer
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public void bufferHeader(int priority, Appendable buffer) throws IOException {
    bufferHeader(mType, mName, mDescription, mColor, priority, buffer);
  }

  /**
   * Buffer header.
   *
   * @param type        the type
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @param priority    the priority
   * @param buffer      the buffer
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void bufferHeader(String type, String name, String description, Color color, int priority,
      Appendable buffer) throws IOException {
    buffer.append(TRACK_PREFIX);
    buffer.append(" type=").append(type);
    buffer.append(" name=\"").append(TextUtils.truncate(name, 15)).append("\"");
    buffer.append(" description=").append(TextUtils.quote(description));
    buffer.append(" priority=").append(Integer.toString(priority));
    UCSCTrack.formatColor(color, buffer);
  }

  /**
   * Format color.
   *
   * @param color  the color
   * @param buffer the buffer
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void formatColor(Color color, Appendable buffer) throws IOException {
    if (color == null) {
      return;
    }

    buffer.append(" color=");
    buffer.append(Integer.toString(color.getRed()));
    buffer.append(TextUtils.COMMA_DELIMITER);
    buffer.append(Integer.toString(color.getGreen()));
    buffer.append(TextUtils.COMMA_DELIMITER);
    buffer.append(Integer.toString(color.getBlue()));
  }

  /**
   * Gets the header.
   *
   * @param type        the type
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @return the header
   */
  public static String getHeader(String type, String name, String description, Color color) {
    return getHeader(type, name, description, color, 1);
  }

  /**
   * Gets the header.
   *
   * @param type        the type
   * @param name        the name
   * @param description the description
   * @param color       the color
   * @param priority    the priority
   * @return the header
   */
  public static String getHeader(String type, String name, String description, Color color, int priority) {
    StringBuilder buffer = new StringBuilder();

    try {
      bufferHeader(type, name, description, color, priority, buffer);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return buffer.toString();
  }

  /**
   * Checks for header.
   *
   * @param file the file
   * @return true, if successful
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static boolean hasHeader(Path file) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    boolean header = false;

    try {
      header = reader.readLine().startsWith(TRACK_PREFIX);
    } finally {
      reader.close();
    }

    return header;
  }

  /**
   * Checks if is track line.
   *
   * @param line the line
   * @return true, if is track line
   */
  public static boolean isTrackLine(String line) {
    return line.startsWith(TRACK_PREFIX);
  }

  /**
   * Gets the name.
   *
   * @param file the file
   * @return the name
   */
  public static String getName(Path file) {
    String name = file.getFileName().toString();

    name = name.substring(0, name.lastIndexOf(TextUtils.PERIOD));

    return name;
  }

  /**
   * Parses the color.
   *
   * @param matcher the matcher
   * @return the color
   */
  public static Color parseColor(Matcher matcher) {
    Color color = new Color(Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(2)),
        Integer.parseInt(matcher.group(3)));

    return color;
  }

  /**
   * Write.
   *
   * @param track the track
   * @param file  the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(UCSCTrack track, Path file) throws IOException {
    write(CollectionUtils.asList(track), file);
  }

  /**
   * Write.
   *
   * @param tracks the tracks
   * @param file   the file
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void write(List<UCSCTrack> tracks, Path file) throws IOException {
    BufferedWriter writer = FileUtils.newBufferedWriter(file);

    int priority = 1;

    try {
      for (UCSCTrack track : tracks) {
        writer.write(track.getHeader(priority));
        writer.newLine();

        for (Entry<Chromosome, List<GenomicElement>> item : track) {
          for (GenomicElement e : item.getValue()) {
            writer.write(e.toString());
            writer.newLine();
          }
        }

        ++priority;
      }
    } finally {
      writer.close();
    }
  }

  /**
   * Gets the track line.
   *
   * @param file the file
   * @return the track line
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static String getTrackLine(Path file) throws IOException {
    BufferedReader reader;

    if (PathUtils.getFileExt(file).equals("gz")) {
      reader = Resources.getGzipReader(file);
    } else {
      reader = FileUtils.newBufferedReader(file);
    }

    String track = null;

    try {
      track = reader.readLine();
    } finally {
      reader.close();
    }

    if (track.startsWith("track")) {
      return track;
    } else {
      return null;
    }
  }

  /**
   * Returns a map of the track attributes.
   *
   * @param file the file
   * @return the track attributes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Map<String, String> getTrackAttributes(Path file) throws IOException {
    return getTrackAttributes(getTrackLine(file));

  }

  /**
   * Gets the track attributes.
   *
   * @param line the line
   * @return the track attributes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Map<String, String> getTrackAttributes(String line) throws IOException {
    if (line == null) {
      return null;
    }

    Map<String, String> map = new TreeMap<String, String>();

    Matcher matcher = ATTRIBUTE_PATTERN.matcher(line);

    while (matcher.find()) {
      map.put(matcher.group(1), TextUtils.removeQuotes(matcher.group(2)));
    }

    return map;
  }

  /**
   * Gets the name from track.
   *
   * @param file the file
   * @return the name from track
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static String getNameFromTrack(Path file) throws IOException {
    return getNameFromTrack(getTrackLine(file));
  }

  /**
   * Creates the bed from track line.
   *
   * @param line the line
   * @return the bed
   */
  public static String getNameFromTrack(String line) {
    if (line == null) {
      return null;
    }

    String name = null;

    Matcher matcher = NAME_PATTERN.matcher(line);

    if (matcher.find()) {
      name = matcher.group(1);
    }

    return name;
  }

  /**
   * Gets the description from track.
   *
   * @param file the file
   * @return the description from track
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static String getDescriptionFromTrack(Path file) throws IOException {
    return getDescriptionFromTrack(getTrackLine(file));
  }

  /**
   * Gets the description from track.
   *
   * @param line the line
   * @return the description from track
   */
  public static String getDescriptionFromTrack(String line) {
    if (line == null) {
      return null;
    }

    String description = null;

    Matcher matcher = DESCRIPTION_PATTERN.matcher(line);

    if (matcher.find()) {
      description = matcher.group(1);
    }

    return description;
  }

  /**
   * Gets the color from track.
   *
   * @param file the file
   * @return the color from track
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Color getColorFromTrack(Path file) throws IOException {
    return getColorFromTrack(getTrackLine(file));
  }

  /**
   * Gets the color from track.
   *
   * @param line the line
   * @return the color from track
   */
  public static Color getColorFromTrack(String line) {
    if (line == null) {
      return null;
    }

    Color color = null;

    Matcher matcher = COLOR_PATTERN.matcher(line);

    if (matcher.find()) {
      color = parseColor(matcher);
    }

    return color;
  }

}
