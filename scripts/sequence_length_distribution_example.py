import gzip

import bionumpy as bnp
import numpy as np
import sys

import bionumpy.io.fastq_buffer
from bionumpy.io.parser import NumpyFileReader


@bnp.streamable()
def get_sequence_lengths(chunk):
    return chunk.sequence.shape[1]


def get_sorted_sequence_length_distribution(file_name):
    f = bnp.open(file_name)
    chunks = f.read_chunks()
    counts = bnp.bincount(get_sequence_lengths(chunks))
    sizes = np.nonzero(counts)[0]
    counts = counts[sizes]
    return sizes, counts


def test():
    get_sorted_sequence_length_distribution("example_data/big.fq.gz")

def _test_2():
    get_sorted_sequence_length_distribution("benchmarks/results/dna_sequences/ENCFF689IPX.fq.gz")

def _test_3():
    filename = "benchmarks/results/dna_sequences/ENCFF689IPX.fq.gz"
    reader = NumpyFileReader(gzip.open(filename), buffer_type=bionumpy.io.fastq_buffer.FastQBuffer)
    reader.set_prepend_mode()
    for chunk in reader.read_chunks():
        pass

if __name__ == "__main__":
    file_name = sys.argv[1]
    out_file = sys.argv[2]

    sizes, counts = get_sorted_sequence_length_distribution(file_name)
    with open(out_file, "w") as f:
        f.write("\n".join(["%d %d" % (size, count) for size, count in zip(sizes, counts)]) + "\n")