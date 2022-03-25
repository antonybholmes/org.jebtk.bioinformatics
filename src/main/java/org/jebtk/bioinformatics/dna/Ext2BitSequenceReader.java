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
package org.jebtk.bioinformatics.dna;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.RepeatMaskType;
import org.jebtk.bioinformatics.genomic.Sequence;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.core.io.FileUtils;

/**
 * Encodes DNA in a 2 bit file representing ACGT. All other characters such as N
 * map to A. Bases are encoded in two bits, so 4 bases per byte. A = 0, C = 1, G
 * = 2, T = 3. Files can be accompanied by a corresponding n
 * 
 *
 * @author Antony Holmes
 *
 */
public class Ext2BitSequenceReader extends ChrSequenceReader {

  /** The m N file map. */
  protected Map<Chromosome, Path> mNFileMap = new HashMap<Chromosome, Path>();

  /** The m mask file map. */
  protected Map<Chromosome, Path> mMaskFileMap = new HashMap<Chromosome, Path>();

  /**
   * Directory containing genome Paths which must be of the form chr.n.txt. Each
   * Path must contain exactly one line consisting of the entire chromosome.
   *
   * @param directory the directory
   */
  public Ext2BitSequenceReader(Path directory) {
    super(directory);
  }

  @Override
  public String getName() {
    return "2bit-ext";
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * edu.columbia.rdf.lib.bioinformatics.genome.GenomeAssembly#getSequence(edu.
   * columbia.rdf.lib.bioinformatics.genome.GenomicRegion, boolean,
   * edu.columbia.rdf.lib.bioinformatics.genome.RepeatMaskType)
   */
  @Override
  public final SequenceRegion getSequence(Genome genome, GenomicRegion region, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {
    Chromosome chr = region.getChr();

    if (!mFileMap.containsKey(chr)) {
      Path file = mFile.resolve(chr + ".dna.2bit");

      mFileMap.put(chr, file);

      file = mFile.resolve(chr + ".n.1bit");

      if (FileUtils.exists(file)) {
        mNFileMap.put(chr, file);
      }

      file = mFile.resolve(chr + ".mask.1bit");

      if (FileUtils.exists(file)) {
        mMaskFileMap.put(chr, file);
      }
    }

    return new SequenceRegion(region,
        getSequence2Bit(mFileMap.get(chr), chr, region.getStart(), region.getEnd(), displayUpper, repeatMaskType));
  }

  /**
   * Gets the sequence4 bit.
   *
   * @param file           the Path
   * @param chr            the chr
   * @param start          the start
   * @param end            the end
   * @param displayUpper   the display upper
   * @param repeatMaskType the repeat mask type
   * @return the sequence4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public Sequence getSequence2Bit(Path file, Chromosome chr, int start, int end, boolean displayUpper,
      RepeatMaskType repeatMaskType) throws IOException {

    int s = start - 1;
    int e = end - 1;

    byte[] buf = getBytes2Bit(file, s, e);

    // how many characters to read
    int l = end - start + 1;

    char[] buffer = new char[l];

    // byte mask;
    int v = 0;

    // the offset to start reading from
    int b = s; // % 4;
    int bi = 0; // b / 4;
    int block;

    for (int i = 0; i < l; ++i) {
      block = b % 4;

      // System.err.println("b " + b + " " + buf[bi] + " " + v + " " + bi + "" +
      // chr +
      // " " + start + " " + end);

      switch (block) {
      case 0:
        v = (buf[bi] >> 6);
        break;
      case 1:
        v = (buf[bi] >> 4);
        break;
      case 2:
        v = (buf[bi] >> 2);
        break;
      default:
        v = buf[bi];
        // We are at the end of a byte so the next read must skip to
        // the next byte in the array
        ++bi;
        break;
      }

      // AND with 3 to get the lowest 2 bits
      v &= 3;

      char c = toChar(v); // , repeatMaskType);

      buffer[i] = c;

      ++b;
    }

    buf = getN(chr, start, end);

    if (buf.length > 0) {
      bi = 0;
      b = s;

      for (int i = 0; i < l; ++i) {
        block = b % 8;

        switch (block) {
        case 0:
          v = (buf[bi] >> 7);
          break;
        case 1:
          v = (buf[bi] >> 6);
          break;
        case 2:
          v = (buf[bi] >> 5);
          break;
        case 3:
          v = (buf[bi] >> 4);
          break;
        case 4:
          v = (buf[bi] >> 3);
          break;
        case 5:
          v = (buf[bi] >> 2);
          break;
        case 6:
          v = (buf[bi] >> 1);
          break;
        default:
          v = buf[bi];
          // We are at the end of a byte so the next read must skip to
          // the next byte in the array
          ++bi;
          break;
        }

        v &= 1;

        if (v == 1) {
          buffer[i] = 'N';
        }

        ++b;
      }
    }

    if (repeatMaskType != RepeatMaskType.UPPERCASE) {
      buf = getMask(chr, start, end);

      if (buf.length > 0) {
        if (repeatMaskType == RepeatMaskType.N) {
          // If mask set, change to 'N'

          bi = 0;
          b = s;

          for (int i = 0; i < l; ++i) {
            block = b % 8;

            switch (block) {
            case 0:
              v = (buf[bi] >> 7);
              break;
            case 1:
              v = (buf[bi] >> 6);
              break;
            case 2:
              v = (buf[bi] >> 5);
              break;
            case 3:
              v = (buf[bi] >> 4);
              break;
            case 4:
              v = (buf[bi] >> 3);
              break;
            case 5:
              v = (buf[bi] >> 2);
              break;
            case 6:
              v = (buf[bi] >> 1);
              break;
            default:
              v = buf[bi];
              // We are at the end of a byte so the next read must skip to
              // the next byte in the array
              ++bi;
              break;
            }

            v &= 1;

            if (v == 1) {
              buffer[i] = 'N';
            }

            ++b;
          }
        } else {
          // If mask set, change letter to lowercase
          bi = 0;
          b = s;

          for (int i = 0; i < l; ++i) {
            block = b % 8;

            switch (block) {
            case 0:
              v = (buf[bi] >> 7);
              break;
            case 1:
              v = (buf[bi] >> 6);
              break;
            case 2:
              v = (buf[bi] >> 5);
              break;
            case 3:
              v = (buf[bi] >> 4);
              break;
            case 4:
              v = (buf[bi] >> 3);
              break;
            case 5:
              v = (buf[bi] >> 2);
              break;
            case 6:
              v = (buf[bi] >> 1);
              break;
            default:
              v = buf[bi];
              // We are at the end of a byte so the next read must skip to
              // the next byte in the array
              ++bi;
              break;
            }

            // System.err.println("v " + v + " " + block);

            v &= 1;

            if (v == 1) {
              buffer[i] = toLower(buffer[i]);
            }

            ++b;
          }
        }
      }
    }

    String dna = new String(buffer);

    if (displayUpper) {
      return Sequence.create(dna);
    } else {
      return Sequence.create(dna.toLowerCase());
    }
  }

  /**
   * Returns the number of Ns in a range.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the n
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public byte[] getN(Chromosome chr, int start, int end) throws IOException {
    Path file = mNFileMap.get(chr);

    if (!FileUtils.exists(file)) {
      return EMPTY_BYTES;
    }

    int s = start - 1;
    int e = end - 1;

    return getBytes1Bit(file, s, e);
  }

  /**
   * Returns the repeat mask for a range.
   *
   * @param chr   the chr
   * @param start the start
   * @param end   the end
   * @return the mask
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private byte[] getMask(Chromosome chr, int start, int end) throws IOException {
    Path file = mMaskFileMap.get(chr);

    if (!FileUtils.exists(file)) {
      return EMPTY_BYTES;
    }

    int s = start - 1;
    int e = end - 1;

    return getBytes1Bit(file, s, e);
  }

  /**
   * Gets the bytes4 bit.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return the bytes4 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static byte[] getBytes2Bit(Path file, int start, int end) throws IOException {
    int sb = start / 4;
    int eb = end / 4;

    // System.err.println(sb + " " + eb);

    return getBytes(file, sb, eb);
  }

  /**
   * Gets the bytes 1 bit.
   *
   * @param file  the file
   * @param start the start
   * @param end   the end
   * @return the bytes 1 bit
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static byte[] getBytes1Bit(Path file, int start, int end) throws IOException {
    int sb = start / 8;
    int eb = end / 8;

    // System.err.println(sb + " " + eb);

    return getBytes(file, sb, eb);
  }

  /*
   * public static byte[] getBytes1Bit(Path file, int start, int end) throws
   * IOException { byte[] buf = null;
   * 
   * RandomAccessFile f = FileUtils.newRandomAccess(file);
   * 
   * Map<Integer, byte[]> cacheMap = new HashMap<Integer, byte[]>();
   * 
   * try { f.seek(BLOCK_SIZE_OFFSET_BYTES); int bases = f.readInt(); int blockSize
   * = f.readInt();
   * 
   * // Byte index int sb = start / 8; int eb = end / 8;
   * 
   * int l = end - start + 1;
   * 
   * // The number of bytes we want to return buf = new byte[l];
   * 
   * int bi = 0;
   * 
   * byte[] blockData = new byte[blockSize];
   * 
   * for (int b = sb; b <= eb; ++b) { // Which block this byte is in int bb = b /
   * blockSize;
   * 
   * 
   * 
   * if (cacheMap.containsKey(bb)) { buf[bi] = cacheMap.get(bb)[b % blockSize]; }
   * else { // Data not cached, so we must seek it in the file
   * f.seek(BLOCKS_OFFSET_BYTES + bb * BLOCK_SIZE_BYTES); int offset =
   * f.readInt();
   * 
   * System.err.println(b + " " + bb + " " + offset + " " + (b % blockSize) + " "
   * + blockSize);
   * 
   * if (offset != -1) { f.seek(offset);
   * 
   * f.read(blockData);
   * 
   * buf[bi] = blockData[b % blockSize];
   * 
   * cacheMap.put(bb, blockData); } }
   * 
   * ++bi; } } finally { f.close(); }
   * 
   * if (buf != null) { return buf; } else { return EMPTY_BYTES; } }
   */
}
