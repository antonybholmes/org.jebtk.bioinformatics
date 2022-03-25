/**
 * Copyright 2018 Antony Holmes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.jebtk.bioinformatics.genomic.geb;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomeService;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicElementsDB;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.core.NameGetter;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.json.Json;
import org.jebtk.core.json.JsonParser;
import org.jebtk.core.text.Join;

/**
 * Encode genes in binary binned format.
 *
 * @author Antony Holmes
 */
public class GEBReader extends GenomicElementsDB implements NameGetter {

  private static final long serialVersionUID = 1L;

  public static final int CHECK = 42;
  public static final byte VERSION = 1;

  /**
   * The Constant INT_BYTES represents the bytes used by a 32bit number (8 * 4)
   */
  public static final int INT_BYTES = 4;

  /** The Constant VERSION_OFFSET. */
  public static final int VERSION_BYTE_OFFSET = INT_BYTES;

  public static final int WINDOW_BYTE_OFFSET = VERSION_BYTE_OFFSET + 1;

  /** The Constant WINDOW_BYTE_OFFSET. */
  public static final int GENES_BYTES_OFFSET = VERSION_BYTE_OFFSET + 1;

  public static final int GENOMIC_TYPE_GENE = 1;
  public static final int GENOMIC_TYPE_TRANSCRIPT = 2;
  public static final int GENOMIC_TYPE_EXON = 4;

  public static final int DOUBLE_BYTES = 8;

  public static final int BLOCK_SEPARATOR = 255;

  // 2^16 - 1
  public static final int MAX_CHILDREN = 65535;
  public static final int MAX_TAGS = 255;
  public static final int MAX_VARCHAR_LENGTH = 255;

  private Genome mGenome;

  private DataReader mDataReader;

  private ElementReader mElementReader;

  // private BinReader mBinReader = null;
  private BTreeReader mBTreeReader = null;

  private RadixReader mRadixReader;
  private Path mDir;
  private int mWindow;
  private Chromosome mChr;
  private final String mPrefix;

  /**
   * Instantiates a new GFB genes.
   *
   * @param genome the genome
   * @param window the window
   * @param dir    the dir
   * @throws IOException
   */
  public GEBReader(Path dir, String prefix, Genome genome, int window) throws IOException {

    // mBinReader = new BinReader(dir, genome, window);
    mRadixReader = new RadixReader(dir, prefix, genome, window);
    mDataReader = new DataReader(dir, prefix, genome, window);
    mElementReader = new ElementReader(mDataReader, dir, prefix, genome, window);

    mDir = dir;
    mPrefix = prefix;
    mGenome = genome;
    mWindow = window;
  }

  @Override
  public String getName() {
    return mPrefix;
  }

  public Path getDir() {
    return mDir;
  }

  @Override
  public Iterable<Genome> getGenomes() {
    return CollectionUtils.asList(mGenome);
  }

  @Override
  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
    List<GenomicElement> elements = new ArrayList<GenomicElement>();

    try {
      _find(region, type, elements);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return overlapping(region, elements);
  }

  /**
   * Find genes in the blocks spanning the coordinates. These are the genes most
   * likely to be overlapping the region of interest. A further test is required
   * to test for overlap. This method is designed to narrow down the list of genes
   * 
   * @param chr
   * @param start
   * @param end
   * @param genes
   * @return
   * @throws IOException
   */
  private void _find(GenomicRegion region, GenomicType type, List<GenomicElement> ret) throws IOException {

    if (mChr == null || !region.mChr.equals(mChr)) {
      if (mBTreeReader != null) {
        mBTreeReader.close();
      }

      // mBinReader = new BinReader(mDir, mPrefix, mGenome, chr, mWindow);
      mBTreeReader = new BTreeReader(mDir, mPrefix, mGenome, region.mChr, mWindow);
      mChr = region.mChr;
    }

    // List<Integer> elementAddresses = mBinReader
    // .elementAddresses(chr, start, end);

    List<Integer> elementAddresses = mBTreeReader.elementAddresses(region.mChr, region.mStart, region.mEnd);

    mElementReader.readElements(elementAddresses, type, ret);
  }

  @Override
  public List<GenomicElement> getElements(Genome genome, String search, GenomicType type) {
    return getElements(search, type, false);
  }

  public List<GenomicElement> getElements(String id, GenomicType type, boolean exact) {
    List<Integer> elementAddresses = new ArrayList<Integer>();

    List<GenomicElement> elements = new ArrayList<GenomicElement>();

    try {
      mRadixReader.elementAddresses(id, exact, elementAddresses);
      mElementReader.readElements(elementAddresses, type, elements);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return elements;
  }

  /**
   * Gets the overlapping genes.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @param genes the genes
   * @return the overlapping genes
   */
  private static List<GenomicElement> overlapping(GenomicRegion region, List<GenomicElement> elements) {
    return overlapping(region, elements, 1);
  }

  /**
   * Gets the overlapping genes.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @param genes the genes
   * @param minBp the min bp
   * @return the overlapping genes
   */
  private static List<GenomicElement> overlapping(GenomicRegion region, List<GenomicElement> elements, int minBp) {
    List<GenomicElement> ret = new ArrayList<GenomicElement>(elements.size());

    for (GenomicElement element : elements) {
      GenomicRegion overlap = GenomicRegion.overlap(region, element);

      if (overlap != null && (overlap.getLength() >= minBp)) {
        ret.add(element);
      }
    }

    return ret;
  }

  @Override
  public void add(GenomicElement element) {
    // Do nothing
  }

  public static final Path getFileName(String type, String prefix) {
    return PathUtils.getPath(Join.on('.').values(prefix, type, "geb").toString());
  }

  public static final Path getFileName(String type, String prefix, Chromosome chr) {
    return PathUtils.getPath(Join.on('.').values(prefix, type, chr, "geb").toString());
  }

  public static final Path getIndexFileName(String prefix) {
    return PathUtils.getPath(Join.on('.').values(prefix, "gei").toString());
  }

  public static byte getStrand(Strand strand) {
    if (Strand.isSense(strand)) {
      return 0;
    } else {
      return 1;
    }
  }

  /**
   * Create a reader from a gei file.
   * 
   * @param file
   * @return
   * @throws IOException
   */
  public static GEBReader loadGEI(Path file) throws IOException {
    Path dir = PathUtils.getDir(file);

    Json json = new JsonParser().parse(file);

    // Settings settings = new Settings().loadIniSettings(file);

    Genome genome = GenomeService.getInstance().get(json.get("genome").getString("name"),
        json.get("genome").getString("build"));

    return new GEBReader(dir, json.getString("name"), genome, json.getInt("window"));
  }
}