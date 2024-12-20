import numpy as np

import pytest

import bionumpy as bnp
from bionumpy.streams.memory_mapping import MemMapEncodedRaggedArray as MemMap
from bionumpy.util.testing import assert_encoded_raggedarray_equal


@pytest.fixture
def filename(data_path):
    return data_path / 'big.fq.gz'


@pytest.fixture
def basename(tmp_path, filename):
    return tmp_path / filename.stem


@pytest.fixture
def loader(filename):
    return lambda: (chunk.sequence for chunk in bnp.open(filename))


def test_create(basename, loader):
    all_sequences = np.concatenate(list(loader()))
    mem_mapped_encoded_ragged_array = MemMap.create(loader, basename)
    assert_encoded_raggedarray_equal(all_sequences, mem_mapped_encoded_ragged_array)


def test_load(basename, loader):
    all_sequences = np.concatenate(list(loader()))
    MemMap.create(loader, basename)
    loaded_mem_mapped_encoded_ragged_array = MemMap.load(basename)
    assert_encoded_raggedarray_equal(all_sequences, loaded_mem_mapped_encoded_ragged_array)


@pytest.mark.skip('Needs pooch to download the file')
def test_big_file(tmp_path):
    import pooch
    url = 'https://github.com/bionumpy/bionumpy-example-data/raw/refs/heads/master/big.fq.gz'
    filename = pooch.retrieve(url, known_hash=None, path=tmp_path)
    basename = tmp_path / 'big'
    n_entries = bnp.count_entries(filename)

    loader = lambda: (chunk.sequence for chunk in bnp.open(filename))
    total_size = sum(chunk.size for chunk in loader())
    mem_mapped = MemMap.create(loader, basename)
    n_sequences = len(mem_mapped)
    assert n_sequences == n_entries
    assert total_size == mem_mapped.size
    chunk_size = 32
    for i in range(0, n_sequences, chunk_size):
        chunk = mem_mapped[i:i + chunk_size]
        assert len(chunk) == min(chunk_size, n_sequences - i)
