import time

import numpy as np

import bionumpy as bnp
from bionumpy import Genome
from bionumpy.alignments import alignment_to_interval
from bionumpy.io.indexed_bam import IndexedBamFile


def multi_intervals(bam_filename='../../example_data/ctcf_chr21-22.bam',
                    interval_filename='../../example_data/ctcf.bed.gz',
                    genome_filename='../../example_data/hg38.chrom.sizes'):
    ib = IndexedBamFile(bam_filename)
    intervals = bnp.open(interval_filename).read()
    t = time.time()
    alignments = ib[intervals]
    print('T =', time.time() - t)
    return alignments
    t = time.time()
    my_entries = old_subset(interval_filename, bam_filename, genome_filename)
    print('T =', time.time() - t)
    assert len(alignments) == len(my_entries)

def main():
    import pooch
    bed_filename = pooch.retrieve('https://www.encodeproject.org/files/ENCFF786YUS/@@download/ENCFF786YUS.bed.gz', None)
    genome_filename = pooch.retrieve('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes', None)
    bam_filename = '/home/knut/Downloads/ENCFF321VRF.bam'
    multi_intervals(bam_filename, bed_filename, genome_filename)

def _test():
    multi_intervals()

def old_subset(interval_filename, bam_filename='../../example_data/ctcf_chr21-22.bam',
               genome_filename='../../example_data/hg38.chrom.sizes'):
    genome = Genome.from_file(genome_filename, filter_function=None)
    genome_mask = genome.read_intervals(interval_filename).get_mask()
    chunks = bnp.open(bam_filename).read_chunks()
    all_intervals = []
    for my_entries in chunks:
        mask = genome_mask[alignment_to_interval(my_entries)].any(axis=-1)
        all_intervals.append(my_entries[mask])
    return np.concatenate(all_intervals)


if __name__ == '__main__':
    main()
