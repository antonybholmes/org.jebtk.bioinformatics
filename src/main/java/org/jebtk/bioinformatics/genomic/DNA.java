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

import org.jebtk.core.settings.SettingsService;

/**
 * Functions related to DNA.
 * 
 * @author Antony Holmes
 *
 */
public class DNA {

  // public static final Path DNA_HOME = AppService.MOD_HOME.resolve("dna");

  /**
   * The constant BASE_A_COLOR.
   */
  public static final Color BASE_A_COLOR = SettingsService.getInstance().getColor("bioinformatics.dna.bases.a.color");

  /**
   * The constant BASE_C_COLOR.
   */
  public static final Color BASE_C_COLOR = SettingsService.getInstance().getColor("bioinformatics.dna.bases.c.color");

  /**
   * The constant BASE_G_COLOR.
   */
  public static final Color BASE_G_COLOR = SettingsService.getInstance().getColor("bioinformatics.dna.bases.g.color");

  /**
   * The constant BASE_T_COLOR.
   */
  public static final Color BASE_T_COLOR = SettingsService.getInstance().getColor("bioinformatics.dna.bases.t.color");

  public static final Color BASE_U_COLOR = BASE_T_COLOR;

  /**
   * The constant BASE_N_COLOR.
   */
  public static final Color BASE_N_COLOR = Color.BLACK;

  /**
   * The constant MEGABASE.
   */
  public static final int MEGABASE = 1000000;

  /**
   * The constant KILOBASE.
   */
  public static final double KILOBASE = 1000;

  public static final String N = "N";

  public static final String LN = "n";

  public static final String U = "U";

  public static final String LU = "u";

  public static final String T = "T";
}
