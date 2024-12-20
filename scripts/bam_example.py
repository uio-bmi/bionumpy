import bionumpy as bnp
from npstructures import ragged_slice


def test_bamquality(filename="example_data/test.bam"):
    # Open the aligments file
    alignments = bnp.open(filename).read()

    # Filter out unmapped reads
    alignments = alignments[alignments.flag & 4 == 0]

    # Extract the first cigar operation for each alignments
    start_cigar = alignments.cigar_op[..., 0]
    
    # Get aligments that start with soft-clip
    start_clipped_alignments = alignments[start_cigar == "s"]
    
    # Get the number of softclipped bases
    n_clipped_bases = start_clipped_alignments.cigar_length[..., 0]
    
    # Extract clipped bases
    clipped_bases = ragged_slice(start_clipped_alignments.sequence,
                                 ends=n_clipped_bases)
    
    # Count bases in softclipped regions
    print(bnp.count_encoded(clipped_bases.ravel()))
    
    # Count bases in whole reads
    print(bnp.count_encoded(alignments.sequence.ravel()))


if __name__ == "__main__":
    test_bamquality()
