import pytest
from npstructures.testing import assert_raggedarray_equal
import numpy as np
import bionumpy as bnp
from bionumpy import as_encoded_array
from bionumpy.util.testing import assert_encoded_raggedarray_equal
from bionumpy.io.strops import (int_to_str, ints_to_strings, join, split, str_to_int, str_equal, str_to_float,
                                float_to_strings, str_to_int_with_missing, str_to_float_with_missing)
from bionumpy.io.strops import _str_equal_two_encoded_ragged_arrays


@pytest.fixture()
def ints():
    return [0, 1, 12, 123]


@pytest.fixture()
def decimal_floats():
    return [1.2, 2.11, 3.123]


@pytest.fixture()
def scientific_floats():
    return ["10.e2", "1.32e3", "1.1e-2", "-1.2e12"]


@pytest.fixture()
def decimal_strings():
    return ["1.2", "2.11", "3.123"]


@pytest.fixture()
def floats():
    return ["2.", ".2e1", "-1.2", "2.11", "3.123", "10.e2", "1.32", "1.1e-2", "-1.2e12", "-2", "100", "100e10"]


@pytest.fixture()
def strings():
    return ["0", "1", "12", "123"]


def test_int_to_string(ints, strings):
    for i, s in zip(ints, strings):
        assert str(int_to_str(i)) == s


def test_multiple_int_to_string(ints, strings):
    for i, s in zip(ints_to_strings(ints), strings):
        assert str(i) == s


def test_str_to_int(ints, strings):
    for i, s in zip(ints, str_to_int(strings)):
        assert i == s


def test_str_to_int_with_missing():
    text = as_encoded_array(['1', '', '', '10', ''])
    np.testing.assert_array_equal(str_to_int_with_missing(text), [1, 0, 0, 10, 0])

def test_str_to_float_with_missing():
    text = as_encoded_array(['1', '', '', '10.1', ''])
    np.testing.assert_array_equal(str_to_float_with_missing(text), [1, np.nan, np.nan, 10.1, np.nan])


def test_str_to_float(floats):
    np.testing.assert_array_almost_equal([float(c) for c in floats], str_to_float(floats))


#@pytest.mark.skip("Inaccurate float")
def test_float_to_str(floats):
    _floats = np.array([float(c) for c in floats])
    ra = float_to_strings(_floats)
    np.testing.assert_array_almost_equal([float(row.to_string()) for row in ra],
                                         _floats)


def test_scientific_str_to_float(scientific_floats):
    np.testing.assert_array_almost_equal([float(c) for c in scientific_floats], str_to_float(scientific_floats))


def test_join(strings):
    seqs = as_encoded_array(strings)
    joined = join(seqs)
    true = as_encoded_array("\t".join(strings))
    assert np.all(joined == true)


def test_split(strings):
    joined = as_encoded_array(",".join(strings))
    splitted = split(joined)
    # assert_raggedarray_equal(splitted, as_encoded_array(strings))


def test_str_equal(strings):
    s = as_encoded_array(strings*2)
    mask = str_equal(s, "12")
    np.testing.assert_array_equal(mask, [False, False, True, False]*2)


def test_str_equal_single(strings):
    assert str_equal(strings[2], "12")
    assert not str_equal(strings[3], "12")


def test_str_equal_two_encoded_ragged_arrays():
    r1 = as_encoded_array(["test1", "2", "test123"])
    r2 = as_encoded_array(["test", "2", "test123"])
    r3 = as_encoded_array(["", "", ""])
    assert np.all(str_equal(r1, r2) == [False, True, True])
    assert np.all(str_equal(r1, r3) == [False, False, False])


def test_chromosome_str_equal(data_path):
    bam = bnp.open(data_path / "test.bam").read()
    bam2 = bnp.open(data_path / "test.bam").read()
    assert np.all(bam.chromosome==bam2.chromosome)


if __name__ == "__main__":
    test_split(["1", "12", "123"])
