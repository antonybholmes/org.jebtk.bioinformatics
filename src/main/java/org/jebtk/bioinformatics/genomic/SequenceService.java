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

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.jebtk.core.collections.DefaultHashMap;
import org.jebtk.core.collections.EntryCreator;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.event.ChangeEvent;
import org.jebtk.core.event.ChangeEventProducer;
import org.jebtk.core.event.ChangeListener;
import org.jebtk.core.event.ChangeListeners;
import org.jebtk.core.settings.SettingsService;

/**
 * Service for extracting DNA/RNA from sequences.
 *
 * @author Antony Holmes
 */
public class SequenceService extends SequenceReader
    implements Iterable<Entry<Genome, SequenceReader>>, ChangeEventProducer {

  /**
   * The Class GenomeAssemblyServiceLoader.
   */
  private static class GenomeAssemblyServiceLoader {

    /** The Constant INSTANCE. */
    private static final SequenceService INSTANCE = new SequenceService();
  }

  /**
   * Gets the single instance of GenomeAssemblyService.
   *
   * @return single instance of GenomeAssemblyService
   */
  public static SequenceService getInstance() {
    return GenomeAssemblyServiceLoader.INSTANCE;
  }

  private final List<SequenceReader> mReaders = new UniqueArrayList<>();

  private final IterMap<Genome, SequenceReader> mGenomeMap = new IterTreeMap<>();

  private final ChangeListeners mListeners = new ChangeListeners();

  private boolean mAutoLoad = true;

  private final IterMap<Character, Color> mColorMap = DefaultHashMap
      .create(SettingsService.getInstance().getColor("bioinformatics.dna.bases.n.color"));

  private final IterMap<Character, ChangeListeners> mListenerMap = DefaultHashMap.create(() -> new ChangeListeners());

  private SequenceReader mCurrent;

  /**
   * Instantiates a new cytobands.
   */
  private SequenceService() {
    updateColors();
  }

  /**
   * Get a base color.
   * 
   * @param base
   * @return
   */
  public Color getBaseColor(char base) {
    return mColorMap.get(base);
  }

  /**
   * Set the color for a given base.
   * 
   * @param base
   * @param color
   */
  public void setBaseColor(char base, Color color) {
    mColorMap.put(base, color);

    // SysUtils.err().println("Set base color", base, color);

    fireChanged(base);

    fireChanged();
  }

  /**
   * Gets the base A color.
   *
   * @return the base A color
   */
  public Color getBaseAColor() {
    return getBaseColor('A');
  }

  /**
   * Sets the base A color.
   *
   * @param color the new base A color
   */
  public void setBaseAColor(Color color) {
    SettingsService.getInstance().update("bioinformatics.dna.bases.a.color", color);

    setBaseColor('A', color);
  }

  /**
   * Gets the base C color.
   *
   * @return the base C color
   */
  public Color getBaseCColor() {
    return getBaseColor('C');
  }

  /**
   * Sets the base C color.
   *
   * @param color the new base C color
   */
  public void setBaseCColor(Color color) {
    SettingsService.getInstance().update("bioinformatics.dna.bases.c.color", color);
    setBaseColor('C', color);
  }

  /**
   * Gets the base G color.
   *
   * @return the base G color
   */
  public Color getBaseGColor() {
    return getBaseColor('G');
  }

  /**
   * Sets the base G color.
   *
   * @param color the new base G color
   */
  public void setBaseGColor(Color color) {
    SettingsService.getInstance().update("bioinformatics.dna.bases.g.color", color);
    setBaseColor('G', color);
  }

  /**
   * Gets the base T color.
   *
   * @return the base T color
   */
  public Color getBaseTColor() {
    return getBaseColor('T');
  }

  /**
   * Sets the base T color.
   *
   * @param color the new base T color
   */
  public void setBaseTColor(Color color) {
    SettingsService.getInstance().update("bioinformatics.dna.bases.t.color", color);
    setBaseColor('T', color);
  }

  /**
   * Gets the base N color.
   *
   * @return the base N color
   */
  public Color getBaseNColor() {
    return getBaseColor('N');
  }

  /**
   * Sets the base N color.
   *
   * @param color the new base N color
   */
  public void setBaseNColor(Color color) {
    SettingsService.getInstance().update("bioinformatics.dna.bases.n.color", color);
    setBaseColor('N', color);
  }

  /**
   * Reset the colors back to their defaults.
   */
  public void reset() {
    SettingsService.getInstance().resetToDefault("bioinformatics.dna.bases.a.color");
    SettingsService.getInstance().resetToDefault("bioinformatics.dna.bases.c.color");
    SettingsService.getInstance().resetToDefault("bioinformatics.dna.bases.g.color");
    SettingsService.getInstance().resetToDefault("bioinformatics.dna.bases.t.color");
    SettingsService.getInstance().resetToDefault("bioinformatics.dna.bases.n.color");

    updateColors();
  }

  /**
   * Update colors.
   */
  private void updateColors() {
    mColorMap.put('A', SettingsService.getInstance().getColor("bioinformatics.dna.bases.a.color"));

    mColorMap.put('C', SettingsService.getInstance().getColor("bioinformatics.dna.bases.c.color"));

    mColorMap.put('G', SettingsService.getInstance().getColor("bioinformatics.dna.bases.g.color"));

    mColorMap.put('T', SettingsService.getInstance().getColor("bioinformatics.dna.bases.t.color"));

    mColorMap.put('N', SettingsService.getInstance().getColor("bioinformatics.dna.bases.n.color"));

    fireAllChanged();
  }

  private void fireAllChanged() {
    fireChanged('A');
    fireChanged('C');
    fireChanged('G');
    fireChanged('T');
    fireChanged('N');

    fireChanged();
  }

  /**
   * Fire that a particular base has changed.
   * 
   * @param base
   */
  private void fireChanged(char base) {
    mListenerMap.get(base).fireChanged(new ChangeEvent(this));
  }

  /**
   * Receive events for a specific base.
   * 
   * @param base
   * @param l
   */
  public void addChangeListener(char base, ChangeListener l) {
    mListenerMap.get(base).addChangeListener(l);
  }

  public void add(SequenceReader reader) {
    mReaders.add(reader);

    mCurrent = reader;
  }

  /**
   * Indicate that the genome references have changed so it they may need to be
   * cached again.
   */
  private void autoLoad() {
    if (mAutoLoad) {
      // One assembly object can load multiple genomes potentially.
      for (SequenceReader reader : mReaders) {
        try {
          for (Genome genome : reader.getGenomes()) {
            mGenomeMap.put(genome, reader);
          }
        } catch (IOException e) {
          e.printStackTrace();
        }
      }

      mAutoLoad = false;
    }
  }

  /**
   * Return the latest reader added.
   * 
   * @return
   */
  public SequenceReader getCurrent() {
    return mCurrent;
  }

  public SequenceReader get(Genome genome) {
    autoLoad();

    return mGenomeMap.get(genome);
  }

  @Override
  public String getName() {
    return "sequence-reader-service";
  }

  @Override
  public Iterator<Entry<Genome, SequenceReader>> iterator() {
    autoLoad();

    return mGenomeMap.iterator();
  }

  public void cache() {
    mAutoLoad = true;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.jebtk.bioinformatics.genome.GenomeAssembly#getSequence(java.lang.
   * String, org.jebtk.bioinformatics.genome.GenomicRegion, boolean,
   * org.jebtk.bioinformatics.genome.RepeatMaskType)
   */
  @Override
  public SequenceRegion getSequence(Genome genome, GenomicRegion region, boolean displayUpper,
      RepeatMaskType repeatMaskType) {

    // Iterate over all assemblies until one works.

    SequenceReader a = get(genome);

    if (a == null) {
      return null;
    }

    SequenceRegion ret = null;

    try {
      ret = a.getSequence(genome, region, displayUpper, repeatMaskType);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return ret;
  }

  @Override
  public void addChangeListener(ChangeListener l) {
    mListeners.addChangeListener(l);
  }

  @Override
  public void removeChangeListener(ChangeListener l) {
    mListeners.removeChangeListener(l);
  }

  private void fireChanged() {
    fireChanged(new ChangeEvent(this));
  }

  @Override
  public void fireChanged(ChangeEvent e) {
    mListeners.fireChanged(e);
  }

}