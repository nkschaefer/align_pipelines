#! /usr/bin/env python3
import sys
import os
import argparse
import re
import math
"""
Given a fasta index (from samtools faidx), splits genome up into specified
number of chunks and outputs each chunk as a BED file.

This approach will never split a chromosome, so if you have a well-assembled
genome, you cannot split into a higher number than the number of non-filtered
chromosomes.

If you have a highly fragmented genome, this will split into groups of contigs.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faidx", "-f", help="Reference genome fasta index (.fai)",
            required=True)
    parser.add_argument("--out_dir", "-d", help="Output directory to write BED files",
            required=True)
    parser.add_argument("--num_chunks", "-n", help="Number of chunks (this is a maximum; \
if higher than the number of non-filtered chromosomes, it will be set to that number)",
            type=int, required=True)
    parser.add_argument("--rm_seqs", "-r", help="Patterns to exclude from sequence names, \
pipe separated. Default = chrUn|random|alt", default="chrUn|random|alt", required=False)
    parser.add_argument("--minlen", "-m", help="Minimum length for a sequence to be included",
            type=int, required=False, default=1000)
    return parser.parse_args()

def write_chunk(curchunk, block_idx, out_dir, faibase):
    f = open('{}/{}.{}.bed'.format(out_dir, faibase, block_idx), 'w')
    for chrom, chromlen in curchunk:
        print("{}\t0\t{}".format(chrom, chromlen), file=f)
    f.close()

def main(args):
    options = parse_args()
    
    chromfilt = re.compile(options.rm_seqs)
    
    totsize = 0
    chroms = {}
    nchroms = 0

    f = open(options.faidx, 'r')
    for line in f:
        line = line.rstrip()
        if len(line) > 0:
            dat = line.split('\t')
            if len(dat) > 2 and not chromfilt.match(dat[0]):
                if int(dat[1]) > options.minlen:
                    chroms[dat[0]] = int(dat[1])
                    totsize += int(dat[1])
                    nchroms += 1
    f.close()
    
    nchunk = options.num_chunks
    if nchroms < options.num_chunks:
        nchunk = nchroms
    

    faibase = '.'.join(options.faidx.split('/')[-1].split('.')[0:-1])

    blocksize = math.ceil(totsize/nchunk)
    block_idx = 1
    blocksum = 0
    curchunk = []
    for chrom, size in sorted(chroms.items(), key=lambda x: -x[1]):
        if blocksum + size >= blocksize:
            # Write out prev block
            if len(curchunk) > 0:
                write_chunk(curchunk, block_idx, options.out_dir, faibase) 
                curchunk = []
                blocksum = 0
                block_idx += 1
        curchunk.append((chrom, size))
        blocksum += size
    if len(curchunk) > 0:
        write_chunk(curchunk, block_idx, options.out_dir, faibase)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
