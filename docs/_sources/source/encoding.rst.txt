.. _encodings:

Encodings
~~~~~~~~~
A central concept in BioNumPy is the *encoding* of data, such as DNA sequence, base qualities, kmers, etc, to memory-efficient data types that are used internally by BioNumPy.

For instance, when asking BioNumPy to store the DNA-sequence `s = bnp.as_encoded_array("ACGT", bnp.DNAEncoding)`, BioNumPy does not store the letters A, C, T and G, but instead uses an efficient numeric representation of them. However, the user does not need to know how this works internally, and will only need to think about the stored sequence as letters. This is why things like `sequence == "A"` works, even though sequence internally is a numeric array.


Encoding data using a specific encoding
---------------------------------------
When reading data with `bnp.open`, BioNumPy automatically encodes your data with a suitable encoding (determined by the file format). If you for some reason have data from other sources, e.g. as strings or list of strings, you may use `bnp.as_encoded_array` to encode your data. In this case, you should make sure you use an encoding suitable for your data (see supported encodings below).

    >>> import bionumpy as bnp
    >>> sequences = ["ACCT", "AcaATA", "ca"]
    >>> bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    encoded_ragged_array(['ACCT',
                          'ACAATA',
                          'CA'], AlphabetEncoding('ACGT'))

Supported encodings
---------------------
These are the most common encodings that can be used:

* `bnp.BaseEncoding`: Can encode any string
* `bnp.DNAEncoding`: Supports A, C, T and G (not N)
* `bnp.encodings.alphabet_encoding.ACTGnEncoding`: Supports N, A, C, T, G
* `bnp.encodings.alphabet_encoding.ACUGEncoding`: Supports A, C, U, G
* `bnp.encodings.alphabet_encoding.RNAEncoding`: Supports A, C, U, G
* `bnp.encodings.alphabet_encoding.AminoAcidEncoding`: Supports all valid AminoAcids


When having already encoded data (advanced usage)
--------------------------------------------------
In some cases, you may have already encoded data that you want to use with BioNumPy. In this case, you can wrap your data in BioNumPy's `EncodedArray` or `EncodedRaggedArray` class, but you will need to be sure that your data is encoded correctly as BioNumPy does not verify this.

For example, if you have already encoded DNA-sequences so that A is 0, C is 1, G is 2 and T is 3, your data is compatible with `bnp.DNAEncoding`:

	>>> import numpy as np
	>>> already_encoded = np.array([0, 1, 2, 3])
	>>> bnp.EncodedArray(already_encoded, bnp.DNAEncoding)
	encoded_array('ACGT', AlphabetEncoding('ACGT'))




