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
package org.jebtk.bioinformatics.pathway;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.bioinformatics.genomic.GeneSymbol;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.TextUtils;

/**
 * Convert between gene ids/symbols.
 */
public class Conversion {

  /**
   * The member map.
   */
  private Map<String, String> mMap;

  /**
   * Cope with alternative names.
   */
  private Map<String, Set<String>> mAltMap;

  /**
   * Cope with synonyms.
   */
  private Map<String, Set<String>> mSynMap;

  /**
   * The member ref seq map.
   */
  private Map<String, Set<String>> mRefSeqMap;

  /**
   * The member entrez map.
   */
  private Map<String, Set<String>> mEntrezMap;

  /**
   * The member symbol map.
   */
  private Map<String, Set<String>> mSymbolMap;

  /**
   * The member ensembl gene map.
   */
  private Map<String, Set<String>> mEnsemblGeneMap;

  /**
   * The member ensembl transcript map.
   */
  private Map<String, Set<String>> mEnsemblTranscriptMap;

  /**
   * The member chr map.
   */
  private HashMap<String, String> mChrMap;

  /**
   * The member strand map.
   */
  private HashMap<String, Character> mStrandMap;

  /**
   * The member start map.
   */
  private HashMap<String, List<Integer>> mStartMap;

  /**
   * The member end map.
   */
  private HashMap<String, List<Integer>> mEndMap;

  /**
   * Create a conversion tool. Gene symbols
   *
   * @param reader the reader
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public Conversion(BufferedReader reader) throws IOException {

    mMap = new HashMap<String, String>();
    mAltMap = new HashMap<String, Set<String>>();
    mSynMap = new HashMap<String, Set<String>>();
    mRefSeqMap = new HashMap<String, Set<String>>();
    mEntrezMap = new HashMap<String, Set<String>>();
    mSymbolMap = new HashMap<String, Set<String>>();
    mEnsemblGeneMap = new HashMap<String, Set<String>>();
    mEnsemblTranscriptMap = new HashMap<String, Set<String>>();
    mChrMap = new HashMap<String, String>();
    mStrandMap = new HashMap<String, Character>();
    mStartMap = new HashMap<String, List<Integer>>();
    mEndMap = new HashMap<String, List<Integer>>();

    try {
      reader.readLine();

      String line;

      while ((line = reader.readLine()) != null) {
        if (Io.isEmptyLine(line)) {
          continue;
        }

        List<String> tokens = TextUtils.tabSplit(line);

        String globalId = tokens.get(0);
        String globalIdBase = globalId.split("\\.")[0];
        String refseq = tokens.get(1);
        String entrez = tokens.get(2);
        String ensemblTranscript = tokens.get(3);
        String ensemblGene = tokens.get(4);
        String symbol = tokens.get(5);
        List<String> previousSymbols = TextUtils.scSplit(tokens.get(6));
        List<String> synonymousSymbols = TextUtils.scSplit(tokens.get(7));
        String chr = tokens.get(8);
        char strand = tokens.get(9).charAt(0);

        // UCSC notation
        int start = Integer.parseInt(tokens.get(10)) + 1;
        int end = Integer.parseInt(tokens.get(11));

        mMap.put(globalId, globalIdBase);

        if (!refseq.equals(TextUtils.NA)) {
          addId(refseq, globalId, mRefSeqMap, chr, strand, start, end);
          addId(refseq, globalIdBase, mRefSeqMap, chr, strand, start, end);
        }

        if (!entrez.equals(TextUtils.NA)) {
          addId(entrez, globalIdBase, mEntrezMap, chr, strand, start, end);
        }

        if (!ensemblGene.equals(TextUtils.NA)) {
          addId(ensemblGene, globalIdBase, mEnsemblGeneMap, chr, strand, start, end);
        }

        if (!ensemblTranscript.equals(TextUtils.NA)) {
          addId(ensemblTranscript, globalId, mEnsemblTranscriptMap, chr, strand, start, end);
          addId(ensemblTranscript, globalIdBase, mEnsemblTranscriptMap, chr, strand, start, end);
        }

        // We can control which symbols are loaded

        if (!symbol.equals(TextUtils.NA)) {
          addId(symbol, globalIdBase, mSymbolMap, chr, strand, start, end);

          for (String ps : previousSymbols) {
            addAltId(ps, globalIdBase);
          }

          for (String ss : synonymousSymbols) {
            addSynId(ss, globalIdBase);
          }
        }
      }
    } finally {
      reader.close();
    }
  }

  /*
   * private void addId(String id, String rdf, Map<String, String> idMap, String
   * chr, char strand, int start, int end) { mMap.put(id, rdf);
   * mMap.put(id.toUpperCase(), rdf);
   * 
   * idMap.put(rdf, id);
   * 
   * mChrMap.put(rdf, chr); mStrandMap.put(rdf, strand);
   * 
   * if (!mStartMap.containsKey(rdf)) { mStartMap.put(rdf, new
   * ArrayList<Integer>()); }
   * 
   * mStartMap.get(rdf).add(start);
   * 
   * if (!mEndMap.containsKey(rdf)) { mEndMap.put(rdf, new ArrayList<Integer>());
   * }
   * 
   * mEndMap.get(rdf).add(end); }
   */

  /**
   * Adds the id.
   *
   * @param id       the id
   * @param globalId the rdf
   * @param idMap    the id map
   * @param chr      the chr
   * @param strand   the strand
   * @param start    the start
   * @param end      the end
   */
  private void addId(String id, String globalId, Map<String, Set<String>> idMap, String chr, char strand, int start,
      int end) {
    mMap.put(id, globalId);
    mMap.put(id.toUpperCase(), globalId);

    if (!idMap.containsKey(globalId)) {
      idMap.put(globalId, new HashSet<String>());
    }

    idMap.get(globalId).add(id);

    mChrMap.put(globalId, chr);
    mStrandMap.put(globalId, strand);

    if (!mStartMap.containsKey(globalId)) {
      mStartMap.put(globalId, new ArrayList<Integer>());
    }

    mStartMap.get(globalId).add(start);

    if (!mEndMap.containsKey(globalId)) {
      mEndMap.put(globalId, new ArrayList<Integer>());
    }

    mEndMap.get(globalId).add(end);
  }

  /**
   * Allow old gene symbols to be mapped to current annotation.
   *
   * @param id  the id
   * @param rdf the rdf
   */
  private void addAltId(String id, String rdf) {

    if (!mAltMap.containsKey(id)) {
      mAltMap.put(id, new HashSet<String>());
    }

    mAltMap.get(id).add(rdf);

    if (!mAltMap.containsKey(id.toUpperCase())) {
      mAltMap.put(id.toUpperCase(), new HashSet<String>());
    }

    mAltMap.get(id.toUpperCase()).add(rdf);
  }

  /**
   * Adds the syn id.
   *
   * @param id  the id
   * @param rdf the rdf
   */
  private void addSynId(String id, String rdf) {

    if (!mSynMap.containsKey(id)) {
      mSynMap.put(id, new HashSet<String>());
    }

    mSynMap.get(id).add(rdf);

    if (!mSynMap.containsKey(id.toUpperCase())) {
      mSynMap.put(id.toUpperCase(), new HashSet<String>());
    }

    mSynMap.get(id.toUpperCase()).add(rdf);
  }

  /**
   * Gets the global ids ids.
   *
   * @param ids           the ids
   * @param splitOnHyphen the split on hyphen
   * @param caseSensitive the case sensitive
   * @param keepMapped    the keep mapped
   * @param keepUnmapped  the keep unmapped
   * @param geneMap       the rdf map
   * @return the rdf ids
   */
  public void getIntermediateGenes(List<String> ids, boolean splitOnHyphen, boolean caseSensitive, boolean keepMapped,
      boolean keepUnmapped, Map<String, Set<IntermediateGene>> geneMap) {
    // Split dashed versions of genes

    // List<String> newIds = new ArrayList<String>();

    for (String id : ids) {
      boolean split = splitOnHyphen;

      List<String> parts = TextUtils.fastSplit(id, TextUtils.DASH_DELIMITER);

      // all parts must be valid gene symbols
      if (split) {
        int splitCount = 0;

        for (String part : parts) {
          if (lookupGlobalId(part).size() > 0) {
            ++splitCount;
          }
        }

        split = splitCount == parts.size();
      }

      // Cannot be split if anti-sense
      if (split) {
        for (String part : parts) {
          if (part.matches("AS\\d+")) {
            split = false;
            break;
          }
        }
      }

      // Cannot be split if name is gene-number
      if (split) {
        for (String part : parts) {
          if (part.matches("^\\d+$")) {
            split = false;
            break;
          }
        }
      }

      if (split && parts.size() > 1) {
        // we can split the gene implying it is a read-through,
        // in which case check that the individual genes
        // overlap the read through

        // Each part other than the first, must be
        // genomically local to the first

        // If the symbol is official this set will contain
        // one item (the official id). If the symbol is alternative
        // the set will contain all possible ids that this
        // symbol could represent

        // The readthrough id
        // String refRDF = mMap.get(id);

        split = false;

        Set<String> globalIds = lookupGlobalId(id);

        for (String globalId : globalIds) {

          String refChr = mChrMap.get(globalId);
          char refStrand = mStrandMap.get(globalId);

          // all parts must be genomically together

          int splitCount = 0;

          for (String part : parts) {
            Set<String> subGlobalIds = lookupGlobalId(part);

            for (String subGlobalId : subGlobalIds) {
              String chr = mChrMap.get(subGlobalId);
              char strand = mStrandMap.get(subGlobalId);

              // System.err.println(refRDF + " " + rdf + " " + chr + " " +
              // strand + " " +
              // refStrand);

              if (chr.equals(refChr) && strand == refStrand) {
                ++splitCount;
                break;
              }
            }
          }

          // If both parts are on the chr and strand, we may
          // be able to split them
          split = splitCount == parts.size();

          if (split) {

            // If we are still good, test that the coordinates
            // of each part overlap the read-through

            List<Integer> refStarts = mStartMap.get(globalId);
            List<Integer> refEnds = mEndMap.get(globalId);

            int overlapCount = 0;

            for (String part : parts) {
              Set<String> subGlobalIds = lookupGlobalId(part);

              boolean found = false;

              for (String subGlobalId : subGlobalIds) {

                List<Integer> starts = mStartMap.get(subGlobalId);
                List<Integer> ends = mEndMap.get(subGlobalId);

                // Test all start and end pairs of the read-through
                // and a part to see if any overlap. If they do, we
                // assume the gene is part of the read-through and
                // can be split
                for (int i = 0; i < refStarts.size(); ++i) {
                  int refStart = refStarts.get(i);
                  int refEnd = refEnds.get(i);

                  for (int j = 0; j < starts.size(); ++j) {
                    int start = starts.get(j);
                    int end = ends.get(j);

                    if ((start >= refStart && start <= refEnd) || (end >= refStart && end <= refEnd)) {
                      found = true;
                      break;
                    }
                  }

                  if (found) {
                    break;
                  }
                }

                if (found) {
                  break;
                }
              }

              if (found) {
                ++overlapCount;
              }
            }

            // Do all parts overlap the read-through?
            split = overlapCount == parts.size();
          }

          // If we have determined it can be split from
          // one combination of ids, there is no point
          // continuing; we can stop.
          if (split) {
            break;
          }
        }
      }

      for (String part : parts) {
        if (part.length() == 0) {
          continue;
        }

        if (!geneMap.containsKey(id)) {
          geneMap.put(id, new HashSet<IntermediateGene>());
        }

        IntermediateGene gene = new IntermediateGene(part);

        String caseId = caseSensitive ? part : part.toUpperCase();

        Set<String> globalIds = lookupGlobalId(caseId);

        if (globalIds.size() > 0) {
          for (String globalId : globalIds) {
            gene.addId(globalId);
          }
        } else {
          if (keepUnmapped) {
            gene.addId(TextUtils.NA); // id);
          }
        }

        geneMap.get(id).add(gene);
      }
    }
  }

  /**
   * Gets the rdf ids.
   *
   * @param id            the id
   * @param splitOnHyphen the split on hyphen
   * @param caseSensitive the case sensitive
   * @param keepMapped    the keep mapped
   * @param keepUnmapped  the keep unmapped
   * @param rdfMap        the rdf map
   * @return the rdf ids
   */
  public void getRdfIds(String id, boolean splitOnHyphen, boolean caseSensitive, boolean keepMapped,
      boolean keepUnmapped, Map<String, Set<IntermediateGene>> rdfMap) {
    getIntermediateGenes(CollectionUtils.asList(id), splitOnHyphen, caseSensitive, keepMapped, keepUnmapped, rdfMap);
  }

  /**
   * Gets the entrez ids.
   *
   * @param rdfIds       the rdf ids
   * @param keepMapped   the keep mapped
   * @param keepUnmapped the keep unmapped
   * @return the entrez ids
   */
  public Map<String, Set<FinalGene>> getEntrezIds(final Map<String, Set<IntermediateGene>> rdfIds, boolean keepMapped,
      boolean keepUnmapped) {
    return getIds(rdfIds, mEntrezMap, keepMapped, keepUnmapped);
  }

  /**
   * Gets the ref seq ids.
   *
   * @param rdfIds       the rdf ids
   * @param keepMapped   the keep mapped
   * @param keepUnmapped the keep unmapped
   * @return the ref seq ids
   */
  public Map<String, Set<FinalGene>> getRefSeqIds(final Map<String, Set<IntermediateGene>> rdfIds, boolean keepMapped,
      boolean keepUnmapped) {
    return getIds(rdfIds, mRefSeqMap, keepMapped, keepUnmapped);
  }

  /**
   * Gets the gene symbols.
   *
   * @param rdfIds       the rdf ids
   * @param keepMapped   the keep mapped
   * @param keepUnmapped the keep unmapped
   * @return the gene symbols
   */
  public Map<String, Set<FinalGene>> getGeneSymbols(final Map<String, Set<IntermediateGene>> rdfIds, boolean keepMapped,
      boolean keepUnmapped) {
    return getIds(rdfIds, mSymbolMap, keepMapped, keepUnmapped);
  }

  /**
   * Gets the ensembl transcripts.
   *
   * @param rdfIds       the rdf ids
   * @param keepMapped   the keep mapped
   * @param keepUnmapped the keep unmapped
   * @return the ensembl transcripts
   */
  public Map<String, Set<FinalGene>> getEnsemblTranscripts(final Map<String, Set<IntermediateGene>> rdfIds,
      boolean keepMapped, boolean keepUnmapped) {
    return getIds(rdfIds, mEnsemblTranscriptMap, keepMapped, keepUnmapped);
  }

  /**
   * Gets the ensembl genes.
   *
   * @param rdfIds       the rdf ids
   * @param keepMapped   the keep mapped
   * @param keepUnmapped the keep unmapped
   * @return the ensembl genes
   */
  public Map<String, Set<FinalGene>> getEnsemblGenes(final Map<String, Set<IntermediateGene>> rdfIds,
      boolean keepMapped, boolean keepUnmapped) {
    return getIds(rdfIds, mEnsemblGeneMap, keepMapped, keepUnmapped);
  }

  /**
   * Gets the ids.
   *
   * @param globalGeneMap the rdf map
   * @param idMap         the id map
   * @param keepMapped    the keep mapped
   * @param keepUnmapped  the keep unmapped
   * @return the ids
   */
  private Map<String, Set<FinalGene>> getIds(final Map<String, Set<IntermediateGene>> globalGeneMap,
      Map<String, Set<String>> idMap, boolean keepMapped, boolean keepUnmapped) {
    // List<FinalGene> newIds = new ArrayList<FinalGene>();

    Map<String, Set<FinalGene>> ret = new HashMap<String, Set<FinalGene>>();

    for (String id : globalGeneMap.keySet()) {
      ret.put(id, new HashSet<FinalGene>());

      for (IntermediateGene gene : globalGeneMap.get(id)) {
        FinalGene newGene = new FinalGene(gene.getName());

        // Most genes map to one rdf id, but in the case of old symbols
        // they may map to more than one

        for (String globalId : gene) {
          boolean found = false;

          if (idMap.containsKey(globalId)) {
            for (String newId : idMap.get(globalId)) {
              newGene.addId(new GlobalGene(globalId, newId));
            }

            found = true;
          } else {
            String base = getBaseId(globalId);

            if (idMap.containsKey(base)) {
              for (String newId : idMap.get(base)) {
                newGene.addId(new GlobalGene(base, newId));
              }

              found = true;
            }
          }

          if (!found && keepUnmapped) {
            newGene.addId(new GlobalGene(globalId, TextUtils.NA));
          }
        }

        ret.get(id).add(newGene);
      }
    }

    return ret;
  }

  /**
   * Gets the locations.
   *
   * @param genes the genes
   * @return the locations
   */
  public List<String> getLocations(List<FinalGene> genes) {
    List<String> locations = new ArrayList<String>();

    for (FinalGene gene : genes) {
      locations.add(TextUtils.scJoin(getLocations(gene)));
    }

    return locations;
  }

  /**
   * Gets the locations.
   *
   * @param gene the gene
   * @return the locations
   */
  public List<String> getLocations(FinalGene gene) {
    List<String> locations = new ArrayList<String>();

    for (GlobalGene g : gene) {
      String id = g.getRdfId();

      // System.err.println("id " + id);

      String location = mChrMap.get(id) + ":" + mStartMap.get(id) + "-" + mEndMap.get(id);

      locations.add(location);
    }

    return locations;
  }

  /**
   * Lookup rdf id.
   *
   * @param id the id
   * @return the sets the
   */
  private Set<String> lookupGlobalId(String id) {
    if (mMap.containsKey(id)) {
      return CollectionUtils.toSet(mMap.get(id));
    } else if (mAltMap.containsKey(id)) {
      return mAltMap.get(id);
    } else if (mSynMap.containsKey(id)) {
      return mSynMap.get(id);
    } else {
      return Collections.emptySet();
    }
  }

  /**
   * Gets the base id.
   *
   * @param id the id
   * @return the base id
   */
  private String getBaseId(String id) {
    return id.split("\\.")[0];
  }

  /**
   * Gets the entrez count.
   *
   * @return the entrez count
   */
  public int getEntrezCount() {
    int count = 0;

    for (String globalId : mEntrezMap.keySet()) {
      count += mEntrezMap.get(globalId).size();
    }

    return count;
  }

  /**
   * Gets the entrez id.
   *
   * @param symbol the symbol
   * @return the entrez id
   */
  public GeneSymbol getEntrezId(String symbol) {
    HashMap<String, Set<IntermediateGene>> geneMap = new HashMap<String, Set<IntermediateGene>>();

    getIntermediateGenes(CollectionUtils.asList(symbol), false, false, true, false, geneMap);

    if (geneMap.size() == 0) {
      return null;
    }

    Map<String, Set<FinalGene>> genes = getRefSeqIds(geneMap, true, false);

    String id = genes.keySet().iterator().next();

    return new GeneSymbol(id, genes.get(id).iterator().next().getName());
  }

}
