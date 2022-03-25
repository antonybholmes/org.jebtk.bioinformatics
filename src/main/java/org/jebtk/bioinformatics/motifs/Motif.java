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
package org.jebtk.bioinformatics.motifs;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jebtk.bioinformatics.BaseCounts;
import org.jebtk.bioinformatics.annotation.Species;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.FormattedTxt;
import org.jebtk.core.text.Splitter;
import org.jebtk.core.text.TextUtils;

/**
 * Represents a motif sequence. Each position consists of the counts for each
 * base (a,c,g,t).
 * 
 * @author Antony Holmes
 *
 */
public class Motif implements Comparable<Motif>, Iterable<BaseCounts>, FormattedTxt {
  // private static final Pattern MOTIF_HEADER_PATTERN =
  // Pattern.compile(">([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)");

  /**
   * The constant JASPAR_HEADER_PATTERN.
   */
  private static final Pattern JASPAR_HEADER_PATTERN = Pattern.compile(">([^\\s]+)\\s(.+)");

  /**
   * The constant USER_DB.
   */
  private static final String USER_DB = "User";

  /**
   * Determine how significant a base must be to be considered part of a triplet.
   */
  private static final double TRIPLET_MIN = 0;

  /**
   * The member name.
   */
  private String mName;

  /**
   * The member id.
   */
  private String mId;

  /**
   * The member gene.
   */
  private String mGene;

  /**
   * The member bg pwm.
   */
  private double mBgPwm;

  /**
   * The member bases.
   */
  private List<BaseCounts> mBases;

  /**
   * The member database.
   */
  private String mDatabase;

  /** The m pwm. */
  private double[][] mPwm;

  /** The m triplets. */
  private ArrayList<Integer> mTriplets;

  /** The m organisms. */
  private List<Species> mOrganisms = new ArrayList<Species>();

  /**
   * Construct a new motif with the given name.
   *
   * @param name  the name
   * @param bases the bases
   */
  public Motif(String name, Collection<BaseCounts> bases) {
    this(name, name, bases);
  }

  /**
   * Instantiates a new motif.
   *
   * @param id    the id
   * @param name  the name
   * @param bases the bases
   */
  public Motif(String id, String name, Collection<BaseCounts> bases) {
    this(id, name, name, bases);
  }

  /**
   * Instantiates a new motif.
   *
   * @param id    the id
   * @param name  the name
   * @param gene  the gene
   * @param bases the bases
   */
  public Motif(String id, String name, String gene, Collection<BaseCounts> bases) {
    this(id, name, gene, null, bases);
  }

  /**
   * Instantiates a new motif.
   *
   * @param id       the id
   * @param name     the name
   * @param gene     the gene
   * @param database the database
   * @param bases    the bases
   */
  public Motif(String id, String name, String gene, String database, Collection<BaseCounts> bases) {
    mId = id;
    mName = name;
    mGene = gene;
    mDatabase = database;

    mBases = new ArrayList<BaseCounts>(bases);

    mBgPwm = Math.pow(0.25, mBases.size());

    mPwm = new double[5][bases.size()];

    for (int i = 0; i < bases.size(); ++i) {
      mPwm[0][i] = mBases.get(i).getA();
      mPwm[1][i] = mBases.get(i).getC();
      mPwm[2][i] = mBases.get(i).getG();
      mPwm[3][i] = mBases.get(i).getT();
      mPwm[4][i] = mBases.get(i).getN();
    }
  }

  /**
   * Instantiates a new motif.
   *
   * @param id   the id
   * @param name the name
   * @param gene the gene
   * @param pwm  the pwm
   */
  public Motif(String id, String name, String gene, double[][] pwm) {
    mId = id;
    mName = name;
    mGene = gene;

    mBgPwm = Math.pow(0.25, pwm[0].length);

    mPwm = pwm;
  }

  /**
   * Gets the pwm.
   *
   * @return the pwm
   */
  public double[][] getPwm() {
    return mPwm;
  }

  /**
   * Gets the counts.
   *
   * @param p the p
   * @return the counts
   */
  public BaseCounts getCounts(int p) {
    return mBases.get(p);
  }

  /**
   * Gets the count.
   *
   * @param base the base
   * @param p    the p
   * @return the count
   */
  public double getCount(char base, int p) {
    return getCounts(p).getCount(base);
  }

  /**
   * Gets the id.
   *
   * @return the id
   */
  public String getId() {
    return mId;
  }

  /**
   * Returns the name of the motif.
   *
   * @return the name
   */
  public String getName() {
    return mName;
  }

  /**
   * Gets the gene.
   *
   * @return the gene
   */
  public String getGene() {
    return mGene;
  }

  /**
   * Gets the database.
   *
   * @return the database
   */
  public String getDatabase() {
    return mDatabase;
  }

  /**
   * Returns the number of bases in the motif.
   *
   * @return the base count
   */
  public int getBaseCount() {
    return mPwm[0].length;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Iterable#iterator()
   */
  @Override
  public Iterator<BaseCounts> iterator() {
    return mBases.iterator();
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    StringBuilder buffer = new StringBuilder();

    try {
      formattedTxt(buffer);
    } catch (IOException e) {
      e.printStackTrace();
    }

    return buffer.toString();
  }

  /**
   * Return the background pwm score for a sequence base for this motif.
   *
   * @return the bg pwm
   */
  public double getBgPwm() {
    return mBgPwm;
  }

  /*
   * (non-Javadoc)
   * 
   * @see org.abh.lib.FormattedTxt#formattedTxt(java.lang.Appendable)
   */
  @Override
  public void formattedTxt(Appendable buffer) throws IOException {
    buffer.append(mName).append(TextUtils.NEW_LINE);

    List<String> values = new ArrayList<String>();

    for (BaseCounts c : this) {
      values.add(Double.toString(c.getA()));
    }

    buffer.append(TextUtils.join(values, TextUtils.TAB_DELIMITER)).append(TextUtils.NEW_LINE);

    values = new ArrayList<String>();

    for (BaseCounts c : this) {
      values.add(Double.toString(c.getC()));
    }

    buffer.append(TextUtils.join(values, TextUtils.TAB_DELIMITER)).append(TextUtils.NEW_LINE);

    values = new ArrayList<String>();

    for (BaseCounts c : this) {
      values.add(Double.toString(c.getG()));
    }

    buffer.append(TextUtils.join(values, TextUtils.TAB_DELIMITER)).append(TextUtils.NEW_LINE);

    values = new ArrayList<String>();

    for (BaseCounts c : this) {
      values.add(Double.toString(c.getT()));
    }

    buffer.append(TextUtils.join(values, TextUtils.TAB_DELIMITER)).append(TextUtils.NEW_LINE);
  }

  /**
   * Parse a matrix file from JASPAR.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parseJaspar(Path file) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    Motif motif = null;

    String line;

    List<Motif> motifs = new ArrayList<Motif>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        Matcher matcher = JASPAR_HEADER_PATTERN.matcher(line);

        if (matcher.find()) {
          String id = matcher.group(1);
          String name = matcher.group(2);

          List<String> a = TextUtils.tabSplit(reader.readLine());
          List<String> c = TextUtils.tabSplit(reader.readLine());
          List<String> g = TextUtils.tabSplit(reader.readLine());
          List<String> t = TextUtils.tabSplit(reader.readLine());

          int l = a.size();

          List<BaseCounts> counts = new ArrayList<BaseCounts>();

          for (int i = 0; i < l; ++i) {
            // Convert the values of each column to percentages

            double af = Double.parseDouble(a.get(i));
            double cf = Double.parseDouble(c.get(i));
            double gf = Double.parseDouble(g.get(i));
            double tf = Double.parseDouble(t.get(i));

            counts.add(new BaseCounts(af, cf, gf, tf, true));
          }

          motif = new Motif(id, name, counts);

          motifs.add(motif);
        }
      }
    } finally {
      reader.close();
    }

    return motifs;
  }

  /**
   * Parses the motif.
   *
   * @param file the file
   * @return the motif
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Motif parseMotif(Path file) throws IOException {
    return parseMotif(file, USER_DB);
  }

  /**
   * Return the first motif from a file.
   *
   * @param file     the file
   * @param database the database
   * @return the motif
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Motif parseMotif(Path file, String database) throws IOException {
    return parseMotifs(file, database).get(0);
  }

  /**
   * Parse a matrix file from JASPAR.
   *
   * @param file     the file
   * @param database the database
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parseMotifs(Path file, String database) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    return parseMotifs(reader, database);
  }

  /**
   * Parses the motifs.
   *
   * @param reader   the reader
   * @param database the database
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parseMotifs(BufferedReader reader, String database) throws IOException {
    Motif motif = null;

    String line;

    List<Motif> motifs = new ArrayList<Motif>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        line = line.trim();

        List<String> tokens = Splitter.onTab().text(line.substring(1));

        String id = tokens.get(0);

        String name = id;

        if (tokens.size() > 1) {
          name = tokens.get(1);
        }

        String gene = id;

        if (tokens.size() > 2) {
          gene = tokens.get(2);
        }

        List<String> a = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> c = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> g = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> t = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());

        int l = a.size();

        List<BaseCounts> counts = new ArrayList<BaseCounts>(l);

        // System.err.println(a);

        for (int i = 0; i < l; ++i) {
          // Convert the values of each column to percentages

          double af = Double.parseDouble(a.get(i));
          double cf = Double.parseDouble(c.get(i));
          double gf = Double.parseDouble(g.get(i));
          double tf = Double.parseDouble(t.get(i));

          counts.add(new BaseCounts(af, cf, gf, tf, true));
        }

        motif = new Motif(id, name, gene, database, counts);

        motifs.add(motif);
      }
    } finally {
      reader.close();
    }

    return motifs;
  }

  /**
   * Parses the pwm motif.
   *
   * @param file the file
   * @return the motif
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Motif parsePwmMotif(Path file) throws IOException {
    return parsePwmMotifs(file).get(0);
  }

  /**
   * Parse a matrix file from JASPAR.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parsePwmMotifs(Path file) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    return parsePwmMotifs(reader);
  }

  /**
   * Parses the pwm motifs.
   *
   * @param reader the reader
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parsePwmMotifs(BufferedReader reader) throws IOException {
    Motif motif = null;

    String line;

    List<Motif> motifs = new ArrayList<Motif>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        line = line.trim();

        String id = line.substring(1);

        List<String> a = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> c = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> g = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> t = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());

        int l = a.size();

        List<BaseCounts> counts = new ArrayList<BaseCounts>(l);

        // First char is the nucleotide so ignore it
        for (int i = 0; i < l; ++i) {
          // Convert the values of each column to percentages

          double af = Double.parseDouble(a.get(i));
          double cf = Double.parseDouble(c.get(i));
          double gf = Double.parseDouble(g.get(i));
          double tf = Double.parseDouble(t.get(i));

          // double total = af + cf + gf + tf;

          // af /= total;
          // cf /= total;
          // gf /= total;
          // tf /= total;

          counts.add(new BaseCounts(af, cf, gf, tf, true));
        }

        motif = new Motif(id, counts);

        motifs.add(motif);
      }
    } finally {
      reader.close();
    }

    return motifs;
  }

  /**
   * Parses the pwm motif.
   *
   * @param file the file
   * @return the motif
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static Motif parsePwm2Motif(Path file) throws IOException {
    return parsePwm2Motifs(file).get(0);
  }

  /**
   * Parse a matrix file from JASPAR.
   *
   * @param file the file
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parsePwm2Motifs(Path file) throws IOException {
    BufferedReader reader = FileUtils.newBufferedReader(file);

    return parsePwm2Motifs(reader);
  }

  /**
   * Parses the pwm motifs.
   *
   * @param reader the reader
   * @return the list
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static List<Motif> parsePwm2Motifs(BufferedReader reader) throws IOException {
    Motif motif = null;

    String line;

    List<Motif> motifs = new ArrayList<Motif>();

    try {
      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        line = line.trim();

        String id = line.substring(1);

        List<String> a = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> c = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> g = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> t = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());
        List<String> n = Splitter.onTab().ignoreEmptyStrings().text(reader.readLine());

        int l = a.size();

        List<BaseCounts> counts = new ArrayList<BaseCounts>(l);

        // First char is the nucleotide so ignore it
        for (int i = 1; i < l; ++i) {
          // Convert the values of each column to percentages

          double af = Double.parseDouble(a.get(i));
          double cf = Double.parseDouble(c.get(i));
          double gf = Double.parseDouble(g.get(i));
          double tf = Double.parseDouble(t.get(i));
          double nf = Double.parseDouble(n.get(i));

          // double total = af + cf + gf + tf;

          // af /= total;
          // cf /= total;
          // gf /= total;
          // tf /= total;

          counts.add(new BaseCounts(af, cf, gf, tf, nf, true));
        }

        motif = new Motif(id, counts);

        motifs.add(motif);
      }
    } finally {
      reader.close();
    }

    return motifs;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Comparable#compareTo(java.lang.Object)
   */
  @Override
  public int compareTo(Motif m) {
    return mName.compareTo(m.mName);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Motif)) {
      return false;
    }

    return compareTo((Motif) o) == 0;
  }

  /**
   * Reverse complement.
   *
   * @param motifs the motifs
   * @return the list
   */
  public static List<Motif> reverseComplement(List<Motif> motifs) {
    List<Motif> ret = new ArrayList<Motif>();

    for (Motif motif : motifs) {
      ret.add(reverseComplement(motif));
    }

    return ret;
  }

  /**
   * Produces the reverse complement of a motif.
   *
   * @param motif the motif
   * @return the motif
   */
  public static Motif reverseComplement(Motif motif) {
    int l = motif.getBaseCount();

    List<BaseCounts> counts = new ArrayList<BaseCounts>(l);

    // for (int i = 0; i < l; ++i) {
    for (int i = l - 1; i >= 0; --i) {
      BaseCounts c = motif.getCounts(i);

      counts.add(new BaseCounts(c.getT(), c.getG(), c.getC(), c.getA(), c.getN()));
    }

    Motif ret = new Motif(motif.mId, motif.mName, motif.mGene, motif.mDatabase, counts);

    return ret;
  }

  /**
   * Sanitize.
   *
   * @param text the text
   * @return the string
   */
  public static String sanitize(String text) {
    return text.replaceAll("\\\\", "_").replaceAll("\\/", "_").replaceAll("\\s+", "_").replaceAll("_+", "_");
  }

  /**
   * Gets the triplets.
   *
   * @return the triplets
   */
  public List<Integer> getTriplets() {
    if (mTriplets == null) {
      mTriplets = new ArrayList<Integer>();

      for (int i = 0; i < 4; ++i) {
        if (mPwm[i][0] >= TRIPLET_MIN) {
          int i1 = (i + 1) * 100;

          for (int j = 0; j < 4; ++j) {
            if (mPwm[j][1] >= TRIPLET_MIN) {
              int i2 = (j + 1) * 10;

              for (int k = 0; k < 4; ++k) {
                if (mPwm[k][2] >= TRIPLET_MIN) {
                  int i3 = k + 1;

                  int triplet = i1 + i2 + i3;

                  mTriplets.add(triplet);
                }
              }
            }
          }
        }
      }
    }

    return mTriplets;
  }

  /**
   * Sets the organisms.
   *
   * @param organisms the new organisms
   */
  public void setOrganisms(List<Species> organisms) {
    mOrganisms.addAll(organisms);
  }

}
