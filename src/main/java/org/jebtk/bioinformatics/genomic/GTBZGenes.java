package org.jebtk.bioinformatics.genomic;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.TokenFunction;

/**
 * Only load genes when requested.
 * 
 * @author antony
 *
 */
public class GTBZGenes extends FixedGapGenes {
  private static final long serialVersionUID = 1L;

  private GeneParser mParser = new GTB2Parser();
  private final Path mFile;
  private final Map<String, Chromosome> mGeneMap = new HashMap<String, Chromosome>();

  public GTBZGenes(Path file, Genome genome) {
    super(genome);

    mFile = file;
  }

  public GTBZGenes(Path file, Genome genome, GeneParser parser) {
    this(file, genome);

    mParser = parser;
  }

  private void geneChrMap() throws IOException {

    final ZipFile zipFile = FileUtils.newZipFile(mFile);

    try {
      ZipEntry entry = zipFile.getEntry("gene_chr.map");

      if (entry != null) {
        BufferedReader reader = FileUtils.newBufferedReader(zipFile, entry);

        try {
          FileUtils.tokenize(new TokenFunction() {
            @Override
            public void parse(final List<String> tokens) {
              String name = GenesDB.sanitize(tokens.get(0));
              Chromosome chr = ChromosomeService.getInstance().guessChr(mFile, tokens.get(1));

              mGeneMap.put(name, chr);
            }
          }).skipHeader(true).tokens(reader);
        } finally {
          reader.close();
        }
      }
    } finally {
      zipFile.close();
    }
  }

  private void autoLoad() {
    try {
      mParser.parse(mFile, mGenome, this);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void autoLoad(Chromosome chr) {
    if (!contains(chr)) {
      try {
        mParser.parse(mFile, mGenome, chr, this);
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  private void autoLoad(String name) {
    // Load a mapping of gene ids to chrs if not done so
    if (mGeneMap.size() == 0) {
      try {
        geneChrMap();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

    Chromosome chr = mGeneMap.get(GenesDB.sanitize(name));

    if (chr != null) {
      autoLoad(chr);
    }
  }

  /*
   * @Override public void autoFindMainVariants() { autoLoad();
   * 
   * super.autoFindMainVariants(); }
   */

  @Override
  public List<GenomicElement> getElements(Genome genome, String id, GenomicType type) {
    autoLoad(id);

    return super.getElements(genome, id, type);
  }

  @Override
  public List<GenomicElement> find(Genome genome, GenomicRegion region, GenomicType type, int minBp) {
    autoLoad(region.mChr);

    return super.find(genome, region, type, minBp);
  }

  @Override
  public List<GenomicElement> closest(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    autoLoad(region.mChr);

    return super.closest(genome, region, type, minBp);
  }

  @Override
  public List<GenomicElement> closestByTss(Genome genome, GenomicRegion region, GenomicType type, int minBp)
      throws IOException {
    autoLoad(region.mChr);

    return super.closestByTss(genome, region, type, minBp);
  }

  /*
   * @Override public Iterable<String> getIds(String type) { autoLoad();
   * 
   * return super.getIds(type); }
   */
}
