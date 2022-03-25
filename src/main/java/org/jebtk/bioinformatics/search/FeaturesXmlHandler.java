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
package org.jebtk.bioinformatics.search;

import java.io.FileNotFoundException;

import org.jebtk.bioinformatics.genomic.Genome;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

// TODO: Auto-generated Javadoc
/**
 * The class FeaturesXmlHandler.
 */
public class FeaturesXmlHandler extends DefaultHandler {

  /**
   * The current genome.
   */
  private Genome currentGenome = Genome.HG18;

  /**
   * The group.
   */
  private String group;

  /*
   * (non-Javadoc)
   * 
   * @see org.xml.sax.helpers.DefaultHandler#startElement(java.lang.String,
   * java.lang.String, java.lang.String, org.xml.sax.Attributes)
   */
  public final void startElement(String uri, String localName, String qName, Attributes attributes)
      throws SAXException {

    if (qName.equals("genome")) {
      if (attributes.getValue("name").equals("hg19")) {
        currentGenome = Genome.HG19;
      } else {
        currentGenome = Genome.HG18;
      }
    } else if (qName.equals("group")) {
      group = attributes.getValue("name");
    } else if (qName.equals("feature")) {

      String name = attributes.getValue("name");
      String filePrefix = attributes.getValue("file");
      String description = attributes.getValue("description");

      try {
        FeaturesServer.getInstance().addCommonRegionBinarySearch(currentGenome, group, name, description, filePrefix);
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    } else {
      // do nothing
    }
  }
}
