import pytest
import bionumpy as bnp

@pytest.fixture
def pysam_install():
    try:
        import pysam
    except ImportError:
        pytest.skip()


def test_indexed_bam(pysam_install, data_path):
    from bionumpy.io import open_indexed
    bam_filepath = data_path/'ctcf_chr21-22.bam'
    bed_filepath = data_path/'ctcf.bed.gz'
    bam = open_indexed(bam_filepath)
    bed = bnp.open(bed_filepath).read()
    alignments = bam[bed]
    assert len(alignments) == 12649