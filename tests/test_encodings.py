import numpy as np
import pytest
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal
from npstructures import RaggedArray
import bionumpy.encoded_array
import bionumpy as bnp
from bionumpy.util.testing import assert_encoded_array_equal, assert_encoded_raggedarray_equal
from bionumpy import as_encoded_array
from bionumpy.io import FastQBuffer
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.datatypes import SequenceEntryWithQuality
from bionumpy.encodings import DigitEncoding, QualityEncoding, CigarEncoding, DigitEncodingFactory, DNAEncoding, \
    ACGTnEncoding
from bionumpy.encoded_array import NumericEncoding, OneToOneEncoding, BaseEncoding


# from bionumpy.encoded_array import OneToOneEncoding

@pytest.fixture
def seqs():
    return ['acgt', 'gct']


def test_alphabet_conversion(seqs):
    dna = as_encoded_array(seqs, ACGTnEncoding)
    assert_encoded_raggedarray_equal(as_encoded_array(dna, DNAEncoding), seqs)


def test_change_encoding_on_encoded_array():
    a = bionumpy.encoded_array.as_encoded_array("ACTG", bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bionumpy.encoded_array.change_encoding(a, bnp.encodings.DNAEncoding)

    assert np.all(a.encoding.decode(a) == b.encoding.decode(b))


def test_change_encoding_on_encoded_ragged_array():
    a = bionumpy.encoded_array.as_encoded_array(["ACTG", "AC"], bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bionumpy.encoded_array.change_encoding(a, bnp.encodings.DNAEncoding)

    assert_encoded_array_equal(a.encoding.decode(a.ravel()), b.encoding.decode(b.ravel()))


def test_as_encoded_on_already_encoded():
    a = bionumpy.encoded_array.as_encoded_array(["ACTG"], bnp.encodings.DNAEncoding)


class CustomEncoding(OneToOneEncoding):
    def _encode(self, data):
        return data + 1

    def _decode(self, data):
        return data - 1


class CustomNumericEncoding(NumericEncoding):
    def _decode(self, data):
        return data + 10

    def _encode(self, data):
        return data - 10


@pytest.mark.parametrize("data", ["test", ["test1", "test2"]])
def test_public_encode_decode_string(data):
    custom_encoding = CustomEncoding()
    encoded = custom_encoding.encode(data)
    decoded = custom_encoding.decode(encoded)
    encoded2 = custom_encoding.encode(decoded)
    assert np.all(encoded == encoded2)


@pytest.mark.parametrize("data", [np.array([1, 2, 3]), RaggedArray([[1, 10], [1]])])
def test_encode_numeric(data):
    custom_encoding = CustomNumericEncoding()
    encoded = custom_encoding.encode(data)
    decoded = custom_encoding.decode(encoded)
    encoded2 = custom_encoding.encode(decoded)
    assert np.all(encoded == encoded2)


@pytest.mark.parametrize("data",
                         ["1234",
                          ["1234", "5678"],
                          np.array([1, 2, 3, 4]),
                          RaggedArray([[1, 2, 3], [4]])])
def test_digit_encoding(data):
    encoding = DigitEncoding
    encoded = encoding.encode(data)
    decoded = encoding.decode(encoded)
    encoded2 = encoding.encode(decoded)
    assert_raggedarray_equal(encoded, encoded2)


@pytest.mark.parametrize("data", ["!!@-^", ["!!@@", "!+"]])
def test_base_quality_encoding(data):
    encoding = QualityEncoding
    encoded = encoding.encode(data)
    decoded = encoding.decode(encoded)
    encoded2 = encoding.encode(decoded)
    assert np.all(encoded2 == encoded)


TestDigitEncoding = DigitEncodingFactory("1")


@bnpdataclass
class Entry:
    a: TestDigitEncoding


def test_numeric_entry():
    encoding = TestDigitEncoding
    a = Entry.single_entry("1234")
    assert np.all(encoding.decode(a.a) == 48 + np.array([1, 2, 3, 4]))


@pytest.mark.parametrize("data", ["!#", ["!#"]])
def test_quality_encoding(data):
    data = as_encoded_array(data, BaseEncoding)
    encoded_data = as_encoded_array(data, QualityEncoding)
    assert np.all(encoded_data.ravel() == [0, 2])


def test_encoding_sequence_entry():
    s = SequenceEntryWithQuality(
        name=as_encoded_array(['headerishere'], BaseEncoding),
        sequence=as_encoded_array(['CTTGTTGA'], BaseEncoding),
        quality=as_encoded_array(['!#!#!#!#'], BaseEncoding),
    )
    correct = [0, 2, 0, 2, 0, 2, 0, 2]

    assert type(s.quality) == RaggedArray
    assert np.all(s.quality.ravel() == correct)
    data = s
    print(data.quality, data.quality.shape)
    # data = FastQBuffer.dataclass.stack_with_ragged(s)
    assert_raggedarray_equal(data.quality, [correct])


@pytest.mark.parametrize("data", ["!!@-^", ["!!@@", "!+"]])
def test_cigar_encoding(data):
    encoding = CigarEncoding
    encoded = encoding.encode(data)
    decoded = encoding.decode(encoded)
    encoded2 = encoding.encode(decoded)
    assert np.all(encoded2 == encoded)


class CustomNumericEncoding(NumericEncoding):
    def _encode(self, data):
        return data - 1

    def _decode(self, data):
        return data + 1


def test_custom_numeric_encoding():
    encoding = CustomNumericEncoding()
    string = "!#!#"
    array = bnp.EncodedArray(np.array([33, 35, 33, 35]), BaseEncoding)
    encoded = encoding.encode(string)
    encoded2 = encoding.encode(array)
    assert_raggedarray_equal(encoded, encoded2)


def test_bool_encoding():
    from bionumpy.encodings.bool_encoding import BoolStringEncoding, bool_string
    encoding = bool_string
    data = ["True", "False", "True"]
    encoded = encoding.encode(data)
    decoded = encoding.decode(encoded)
    assert np.all(decoded == data)
    print(as_encoded_array(data, encoding))