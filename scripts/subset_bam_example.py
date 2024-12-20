import time


import bionumpy as bnp
from bionumpy.alignments import alignment_to_interval


def unique_intersect(filename_a: str, filename_b: str, chrom_sizes_file: str, out_filename: str):
    # Load the genome
    genome = bnp.Genome.from_file(chrom_sizes_file, filter_function=None)

    # Load the intervals from the second file and get the binary mask
    genome_mask = genome.read_intervals(filename_b).get_mask()

    with bnp.open(out_filename, 'w') as f:
        for i, chunk in enumerate(bnp.open(filename_a).read_chunks()):
            intervals = alignment_to_interval(chunk)
            mask = genome_mask[intervals].any(axis=-1)
            f.write(chunk[mask])


def big_example():
    import pooch
    bed_filename = pooch.retrieve('https://www.encodeproject.org/files/ENCFF786YUS/@@download/ENCFF786YUS.bed.gz', None)
    bam_filename = pooch.retrieve('https://www.encodeproject.org/files/ENCFF321VRF/@@download/ENCFF321VRF.bam', None)
    chrom_sizes = pooch.retrieve('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes', None)
    out_filename = 'tmp.bam'
    t = time.time()
    unique_intersect(bam_filename, bed_filename, chrom_sizes, out_filename)
    print('T =', time.time()-t)


def test_subset_bam():
    out_filename = 'example_data/intersect.bam'
    unique_intersect('example_data/ctcf_chr21-22.bam', 'example_data/ctcf.bed.gz',
                     'example_data/hg38.chrom.sizes',
                     out_filename)
    assert bnp.count_entries(out_filename) == 12649

if __name__ == '__main__':
    import sys
    unique_intersect(*sys.argv[1:])
