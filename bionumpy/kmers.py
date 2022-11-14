import numpy as np

from .encodings.kmer_encodings import KmerEncoding
from .rollable import RollableFunction
from .encodings import DNAEncoding
from .encodings.alphabet_encoding import AlphabetEncoding
from .encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from npstructures.bitarray import BitArray
from .util import convolution, is_subclass_or_instance
import logging
logger = logging.getLogger(__name__)



class KmerEncoder(RollableFunction):
    def __init__(self, k, alphabet_encoding):
        self.window_size = k
        self._k = k
        self._encoding = alphabet_encoding
        self._alphabet_size = alphabet_encoding.alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: EncodedArray) -> int:
        sequence = as_encoded_array(sequence, target_encoding=self._encoding)
        return EncodedArray(sequence.data.dot(self._convolution), KmerEncoding(self._encoding, self._k))

    def inverse(self, kmer_hash: int) -> EncodedArray:
        return EncodedArray(((kmer_hash[:, np.newaxis] // self._convolution) % self._alphabet_size), self._encoding)

    def sample_domain(self, n):
        return EncodedArray((np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(n, self._k)), self._encoding)


def get_kmers(sequence: EncodedRaggedArray, k: int) -> EncodedArray:
    """
    Get kmers for sequences.
    Sequences should be encoded with an AlphabetEncoding (e.g. DNAEncoding).

    Parameters
    ----------
    sequence : EncodedRaggedArray
        Sequences to get kmers from
    k : int
        The kmer size (1-31)

    Returns
    -------
    EncodedRaggedArray
        Kmers from the sequences.

    Examples
    --------
    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "AAA", "TTGGC"], bnp.DNAEncoding)
    >>> bnp.kmers.get_kmers(sequences, 3)
    encoded_ragged_array([[ACT, CTG],
                          [AAA],
                          [TTG, TGG, GGC]], 3merEncoding(AlphabetEncoding('ACGT')))
    """

    assert 0 < k < 32, "k must be larger than 0 and smaller than 32"
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding), \
        "Sequence needs to be encoded with an AlphabetEncoding, e.g. DNAEncoding"

    result = KmerEncoder(k, sequence.encoding).rolling_window(sequence)
    return result


@convolution
def fast_hash(sequence, k):
    """
    A faster veresion of get_kmers. Works only for data already encoded with
    an AlphabetEncoding.
    """
    assert isinstance(sequence, EncodedArray), sequence
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding)

    bit_array = BitArray.pack(sequence.data, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    output_encoding = KmerEncoding(sequence.encoding, k)
    return EncodedArray(hashes, output_encoding)


