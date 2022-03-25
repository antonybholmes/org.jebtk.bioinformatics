/**
 * Copyright 2017 Antony Holmes
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
package org.jebtk.bioinformatics.genomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.IterHashMap;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.json.Json;
import org.jebtk.core.json.JsonParser;
import org.jebtk.core.text.TextUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The Class ChromosomeParser.
 */
public class ChromosomeDirs extends GenomeDirs implements ChromosomeReader, Iterable<Chromosome> {

  private static final Logger LOG = LoggerFactory.getLogger(ChromosomeDirs.class);

  public static final String EXT = "chrs.gz";

  private Genome mGenome;

  private IterMap<Integer, Chromosome> mChrIdMap = new IterHashMap<Integer, Chromosome>();

  private IterMap<String, Chromosome> mChrMap = new IterHashMap<String, Chromosome>();

  private List<Chromosome> mChrs = new ArrayList<Chromosome>();

  private Map<Chromosome, Integer> mSizeMap = new HashMap<Chromosome, Integer>();

  private boolean mAutoLoad = true;

  public ChromosomeDirs(Genome genome, Collection<Path> dirs) {
    super(dirs);

    mGenome = genome;
  }

  public ChromosomeDirs(Genome genome, Path... dirs) {
    super(dirs);

    mGenome = genome;
  }

  private void autoLoad() throws IOException {
    if (mAutoLoad) {
      mAutoLoad = false;

      for (Path dir : mDirs) {
        // Path genomeDir = dir.resolve(mGenome.getName());

        LOG.info("Looking for data in {} {}", dir, mGenome);

        // if (FileUtils.isDirectory(genomeDir)) {

        // Path dbDir = genomeDir.resolve(mGenome.getBuild());

        // LOG.info("Looking for chromosomes in {}", dbDir, mGenome);

        // if (FileUtils.isDirectory(dbDir)) {
        List<Path> files = FileUtils.endsWith(dir, EXT);

        for (Path file : files) {
          load(file);
        }
      }
      // }
      // }

      Collections.sort(mChrs);

    }
  }

  private void load(Path file) throws IOException {
    LOG.info("Discovered chromosome info in {}.", file);

    load(file, this);
  }

  public void cache() {
    mAutoLoad = true;
  }

  /**
   * Gets the id.
   *
   * @param chr the chr
   * @return the id
   */
  public int getId(String chr) {
    Chromosome c = chr(chr);

    if (c != null) {
      return c.getId();
    } else {
      return -1;
    }
  }

  @Override
  public Iterator<Chromosome> iterator() {
    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return mChrs.iterator();
  }

  @Override
  public Chromosome chr(String chr) {
    // LOG.info("Request {} {}", chr, getMapId(chr));

    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    String fc = Chromosome.getShortName(chr);

    if (!mChrMap.containsKey(fc)) {
      add(Chromosome.newChr(chr));
    }

    Chromosome ret = mChrMap.get(fc);

    // if (ret.getName().contains("_")) {
    // LOG.info("Request _ {} {} {}", chr, fc, ret);
    // }

    // if (ret.getShortName().contains("_")) {
    // LOG.info("Short Request _ {} {} {}", chr, fc, ret);
    // }

    return ret;
  }

  public Chromosome chr(int id) {
    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return mChrIdMap.get(id);
  }

  @Override
  public int size(Chromosome chr) {
    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return mSizeMap.getOrDefault(chr, -1);
  }

  protected void add(Chromosome chr) {
    mChrIdMap.put(chr.getId(), chr);
    // mChrMap.put(Integer.toString(chr.getId()), chr);
    mChrMap.put(chr.getShortName().toUpperCase(), chr);
    mChrs.add(chr);

    // mAutoLoad = true;
  }

  /**
   * Returns the genome reference, for example hg19.
   * 
   * @return
   */
  @Override
  public Genome getGenome() {
    return mGenome;
  }

  public Chromosome randChr() {
    try {
      autoLoad();
    } catch (IOException e) {
      e.printStackTrace();
    }

    List<String> ids = CollectionUtils.toList(mChrMap.keySet());

    return mChrMap.get(ids.get(Mathematics.rand(ids.size())));
  }

  /*
   * public static ChromosomeDirs parse(Path file) throws IOException {
   * LOG.info("Reading chromosome info from {}", file);
   * 
   * BufferedReader reader = FileUtils.newBufferedReader(file);
   * 
   * ChromosomeDirs ret = null;
   * 
   * try { ret = parse(reader); } finally { reader.close(); }
   * 
   * return ret; }
   */

  /*
   * private static ChromosomeDirs parse(BufferedReader reader) throws IOException
   * {
   * 
   * // The first token contains the names etc, ignore the rest of the line String
   * species = TextUtils.tabSplit(reader.readLine()).get(0); String g =
   * TextUtils.tabSplit(reader.readLine()).get(0);
   * 
   * Genome genome = GenomeService.getInstance().genome(species, g);
   * 
   * ChromosomeDirs ret = new ChromosomeDirs(genome);
   * 
   * // Skip header reader.readLine();
   * 
   * String line; List<String> tokens;
   * 
   * while ((line = reader.readLine()) != null) { tokens =
   * TextUtils.tabSplit(line);
   * 
   * // int id = Integer.parseInt(tokens.get(0)); String name = tokens.get(1); int
   * size = Integer.parseInt(tokens.get(2));
   * 
   * Chromosome chr = Chromosome.newChr(name, genome, size);
   * 
   * ret.add(chr); }
   * 
   * return ret; }
   */

  /**
   * Parses the.
   *
   * @param file the file
   * @return the chromosome sizes
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static ChromosomeDirs parseJson(Path file) throws IOException {

    Json json = new JsonParser().parse(file);

    Genome genome = GenomeService.getInstance().guessGenome(json.getString("genome"));

    ChromosomeDirs ret = new ChromosomeDirs(genome);

    Json chrsJson = json.get("chromosomes");

    for (int i = 0; i < chrsJson.size(); ++i) {
      Json chrJson = chrsJson.get(i);

      int id = chrJson.getInt("id");
      String name = chrJson.getString("name");
      int size = chrJson.getInt("size");

      Chromosome chr = Chromosome.newChr(name);

      ret.mChrIdMap.put(id, chr);
      ret.mChrMap.put(Integer.toString(id), chr);
      ret.mChrMap.put(chr.getShortName().toUpperCase(), chr);
      ret.mSizeMap.put(chr, size);
    }

    return ret;
  }

  private static void load(Path file, ChromosomeDirs ret) throws IOException {

    Genome g = ret.getGenome();
    String species = g.getName().toLowerCase();

    LOG.info("Loading chrs from {} {}...", file, g);

    BufferedReader reader = FileUtils.newBufferedReader(file);

    try {
      // Skip header
      reader.readLine();
      reader.readLine();
      reader.readLine();

      String line;
      List<String> tokens;

      while ((line = reader.readLine()) != null) {
        tokens = TextUtils.tabSplit(line);

        // int id = Integer.parseInt(tokens.get(0));
        String name = tokens.get(1);

        if (name.contains("_")) {
          continue;
        }

        int size = Integer.parseInt(tokens.get(2));

        Chromosome chr;

        switch (species) {
        case "human":
          chr = Chromosome.newHumanChr(name);
          break;
        case "mouse":
          chr = Chromosome.newMouseChr(name);
          break;
        default:
          chr = Chromosome.newChr(name);
          break;
        }

        // System.err.println(name + " " + chr);
        ret.mSizeMap.put(chr, size);
        ret.add(chr);
      }
    } finally {
      reader.close();
    }

    LOG.info("Finished.");
  }

}
