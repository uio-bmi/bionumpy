import time

import bionumpy as bnp


def unique_intersect(filename_a: str, filename_b: str, chrom_sizes_file: str, out_filename: str):
    # Load the genome
    genome = bnp.Genome.from_file(chrom_sizes_file, filter_function=None)

    # Load the intervals from the second file and get the binary mask
    genome_mask = genome.read_intervals(filename_b).get_mask()

    with bnp.open(out_filename, 'w') as f:
        for intervals in bnp.open(filename_a).read_chunks():
            mask = genome_mask[intervals].any(axis=-1)
            f.write(intervals[mask])


def test():
    out_filename = 'example_data/ctcf_intersect.bed.gz'
    unique_intersect('example_data/ctcf.bed.gz', 'example_data/znf263.bed.gz',
                     'example_data/hg38.chrom.sizes',
                     out_filename)

    assert bnp.count_entries(out_filename) == 3951

if __name__ == '__main__':
    import sys
    t = time.time()
    unique_intersect(*sys.argv[1:])
    print("Time:", time.time()-t)