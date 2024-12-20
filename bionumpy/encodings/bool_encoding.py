from typing import List

import numpy as np

from bionumpy.encoded_array import Encoding, EncodedArray
from bionumpy.encodings.string_encodings import StringEncoding


class BoolStringEncoding(Encoding):
    '''
    >>> from bionumpy.encodings.bool_encoding import bool_string
    >>> bool_string.encode(['True', 'False', 'True'])
    array([ True, False,  True])
    >>> bool_string.decode([False, False])
    encoded_ragged_array(['False',
                          'False'])
    '''
    returns_raw = True
    def __init__(self, true_string: str = 'True', false_string: str = 'False'):
        self._true_string = true_string
        self._false_string = false_string
        self._string_encoding = StringEncoding([false_string, true_string])
        self._lookup = np.array([false_string, true_string])

    def get_labels(self) -> List[str]:
        return [self._false_string, self._true_string]

    def encode(self, encoded_ragged_array):
        s = self._string_encoding.encode(encoded_ragged_array)
        return s.raw().astype(bool)

    def decode(self, encoded_array):
        a = EncodedArray(np.asanyarray(encoded_array).astype(int), self._string_encoding)
        return self._string_encoding.decode(a)

bool_string = BoolStringEncoding()