import numpy as np


def hash_function(offset, mod):
    def f(kmer):
        return (kmer ^ offset) % mod
    return f


class BloomFilter:
    def __init__(self, hash_functions, mask):
        self._hash_functions = hash_functions
        self._mask = mask

    @classmethod
    def from_hash_functions_and_seqeuences(cls, hash_functions, sequence, mask_size):
        mask = np.zeros(mask_size, dtype=bool)
        for function in hash_functions:
            mask[function(sequence)] = True
        return cls(hash_functions, mask)

    def __getitem__(self, idx):
        return np.all([self._mask[h(idx)] for h in self._hash_functions],
                      axis=0)
