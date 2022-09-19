import pytest
import numpy as np
from io import BytesIO
from bionumpy.files import NumpyFileReader, NpDataclassReader
from .buffers import buffer_texts, combos, big_fastq_text


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta"])
def test_buffer_read(buffer_name):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read()
    assert list(data) == true_data


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta"])
def test_buffer_read_chunks(buffer_name):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = np.concatenate(list(NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read_chunks()))
    assert list(data) == true_data


def test_read_big_fastq():
    io_obj = BytesIO(bytes(big_fastq_text*1000, encoding="utf8"))
    for b in NpDataclassReader(NumpyFileReader(io_obj, combos["fastq"][2], chunk_size=1000)).read_chunks():
        print(b)
