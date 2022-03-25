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

import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.BaseCounts;
import org.jebtk.bioinformatics.annotation.Species;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

// TODO: Auto-generated Javadoc
/**
 * The class KeyXmlHandler.
 */
public class MotifXmlHandler extends DefaultHandler {

  /** The m motifs. */
  private List<Motif> mMotifs = new ArrayList<Motif>();

  /** The m name. */
  private String mName;

  /** The m id. */
  private String mId;

  /** The m organisms. */
  private List<Species> mOrganisms = new ArrayList<Species>();

  /** The m counts. */
  private List<BaseCounts> mCounts = new ArrayList<BaseCounts>();

  /** The m db name. */
  private String mDbName;

  /** The m ret. */
  private Motifs mRet;

  // mPwm = new double[4][bases.size()];

  /*
   * (non-Javadoc)
   * 
   * @see org.xml.sax.helpers.DefaultHandler#startElement(java.lang.String,
   * java.lang.String, java.lang.String, org.xml.sax.Attributes)
   */
  public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {

    if (qName.equals("motifs")) {
      mDbName = attributes.getValue("name");
    } else if (qName.equals("motif")) {
      mName = attributes.getValue("name");
      mId = attributes.getValue("id");
    } else if (qName.equals("pos") || qName.equals("position")) {
      double a = Double.parseDouble(attributes.getValue("score-a"));
      double c = Double.parseDouble(attributes.getValue("score-c"));
      double g = Double.parseDouble(attributes.getValue("score-g"));

      double t;

      if (attributes.getValue("score-t") != null) {
        t = Double.parseDouble(attributes.getValue("score-t"));
      } else if (attributes.getValue("score-u") != null) {
        t = Double.parseDouble(attributes.getValue("score-u"));
      } else {
        t = -1;
      }

      mCounts.add(new BaseCounts(a, c, g, t, true));
    } else {

    }

  }

  /*
   * (non-Javadoc)
   * 
   * @see org.xml.sax.helpers.DefaultHandler#endElement(java.lang.String,
   * java.lang.String, java.lang.String)
   */
  public void endElement(String uri, String localName, String qName) throws SAXException {

    if (qName.equals("motifs")) {
      mRet = new Motifs(mDbName, mMotifs);
    } else if (qName.equals("motif")) {
      Motif motif = new Motif(mId, mName, mName, mDbName, mCounts);

      motif.setOrganisms(mOrganisms);

      mMotifs.add(motif);

      mCounts.clear();
      mOrganisms.clear();
    } else {
      // Do nothing
    }
  }

  /**
   * Gets the motifs.
   *
   * @return the motifs
   */
  public Motifs getMotifs() {
    return mRet;
  }
}
