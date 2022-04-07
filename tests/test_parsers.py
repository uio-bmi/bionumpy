import numpy as np
from bionumpy.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.kmers import TwoBitHash
from bionumpy.bed_parser import Interval, BedBuffer, VCFBuffer, SNP

from .buffers import fastq_buffer, twoline_fasta_buffer, bed_buffer, vcf_buffer

def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)

def test_twoline_fasta_buffer(twoline_fasta_buffer):
    buf = TwoLineFastaBuffer.from_raw_buffer(twoline_fasta_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]

def test_fastq_buffer(fastq_buffer):
    buf = FastQBuffer.from_raw_buffer(fastq_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]

def test_bed_buffer(bed_buffer):
    buf = BedBuffer.from_raw_buffer(bed_buffer)
    intervals = buf.get_intervals()
    assert intervals == Interval(
        [[1, 3],
         [40, 60],
         [400, 600]])

def test_vcf_buffer(vcf_buffer):
    buf = VCFBuffer.from_raw_buffer(vcf_buffer)
    snps = list(buf.get_snps())
    true = [SNP("chr1",	88362,	"A",	"G"),
            SNP("chr1",	887560,	"A",	"C"),
            SNP("chr2",	8878,	"A",	"C")]

