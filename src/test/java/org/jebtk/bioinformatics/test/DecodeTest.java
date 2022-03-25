/**
 * Copyright 2016 Antony Holmes
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
package org.jebtk.bioinformatics.test;

import java.io.IOException;

import org.junit.Test;

public class DecodeTest {
  @Test
  public void posTest() throws IOException {

    /*
     * GFBGenes g = new GFBGenes("gencode", "grch38", 1000,
     * PathUtils.getPath("/home/antony/Desktop/gff/gfb/grch38/"));
     * 
     * List<GenomicEntity> r = g.findGenes("chr3:187721377-187745727");
     * //List<GenomicEntity> r = g.findGenes("chrM:1-10000");
     * 
     * System.err.println("g " + r.size());
     * 
     * for (GenomicEntity gene : r) { System.err.println(gene + " " + gene.getTags()
     * + " " + gene.getChildCount(GenomicType.TRANSCRIPT));
     * 
     * for (GenomicEntity e : gene.getChildren(GenomicType.TRANSCRIPT)) {
     * System.err.println("exon " + e); } }
     */

  }

//	@Test
//	public void searchTest() throws IOException {
//
//		String search = "bcl6";
//
//		GFBGenes genes = new GFBGenes(new Genome("human", "grch38"), 1000,
//				PathUtils.getPath("/home/antony/Desktop/gff/gfb/grch38/"));
//
//		// GeneSearchQuery query = new GeneSearchQuery(genes);
//
//		List<GenomicElement> r = genes.getGenes(search, GenomicType.TRANSCRIPT);
//
//		for (GenomicElement gene : r) {
//			System.err.println("blo " + gene + " " + gene.getType() + " " + gene.getChildTypes());
//
//			for (GenomicElement e : gene.getChildren(GenomicType.EXON)) {
//				System.err.println(e);
//			}
//		}
//	}
}
