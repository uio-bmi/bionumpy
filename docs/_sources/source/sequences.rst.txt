.. _sequences:

=========
Sequences
=========

`EncodedArray`
==============

Sequence data in BioNumPy are are represented by `EncodedArray` objects. These are basically `numpy` arrays of integers, that have an encoding that specifies which character each integer represent. This representataion allows us to do fast `numpy` operation on the sequences, while still allowing for human readable representation of them. The easiest way to create an `EncodedArray` is to use the `as_encoded_array` function.

    >>> import bionumpy as bnp
    >>> encoded_array = bnp.as_encoded_array("actggtcc")
    >>> encoded_array
    encoded_array('actggtcc')

We see that the encoded array represents the text "actggtcc", but under the hood these are just integers with an ASCII encoding. You seldom have to think about this internal representation, but it is what allow us to write for instance:

    >>> print(encoded_array == "g")
    [False False False  True  True False False False]
    
And get numpy-fast performance for the query. We can also use numpy-like indexing on encoded arrays, so that we can for instance trim the first and last two characters from the sequence:

    >>> encoded_array[2:-2]
    encoded_array('tggt')
    
`EncodedRaggedArray`
====================
When working with multiple sequences we usually have to use `EncodedRaggedArray` objects. These are much like `EncodedArray` objects, but instead uses `npstructures.RaggedArray` to store the integers. This allows us to store seqeunces of differing lengths. The easiest way to create an `EncodedRaggedArray` is to use the `as_encoded_array` function on a list of strings:

    >>> encoded_ragged_array = bnp.as_encoded_array(["ctt", "actg", "ag"])
    >>> encoded_ragged_array
    encoded_ragged_array(['ctt',
                          'actg',
                          'ag'])

These objects also behave very much like numpy arrays, in indexing and broadcasting. For instance, to get the 2nd character of the first and third seqeunce:

    >>> encoded_ragged_array[[0, 2], 1]
    encoded_array('tg')

Absolute basics of `EncodedRaggedArray`
=======================================

As shown above, sequences are represented as `EncodedRaggedArray` objects. Below are some absolute basics of `EncodedRaggedArray` objects:

For this example, instead of creating a `EncodedRaggedArray` object directly, we will simulate some sequences using the
BioNumpy's `simulate_sequences` function:

    >>> import bionumpy as bnp
    >>> import numpy as np
    >>> from bionumpy.simulate import simulate_sequences
    >>> rng = np.random.default_rng(seed=1)
    >>> np.random.seed(1)
    >>> aa_alphabet = bnp.encodings.alphabet_encoding.AminoAcidEncoding.get_labels()
    >>> named_seqs = simulate_sequences(aa_alphabet, {f's{i}': np.random.randint(5,20) for i in range(10)}, rng)
    >>> my_seqs = named_seqs.sequence
    >>> named_seqs # print the sequences
     SequenceEntry with 10 entries
                         name                 sequence
                           s0               LMSYAEVYGH
                           s1         WKGVGKQNCAWSVNVH
                           s2        LTDHDL*DKKWFMGASC
                           s3            GMMD*S*CSHNYG
                           s4           SEH*KMHDKQLTIP
                           s5         TYKASNWLICLQTVFP
                           s6               TGIVPMRM*S
                           s7                    CENVC
                           s8                    RSTWF
                           s9                   NTIFMC



Indexing and slicing of `EncodedRaggedArray` objects
-----------------------------------------------------

We can index and slice `EncodedRaggedArray` objects in a similar way to numpy arrays. Below are some examples:

    >>> my_seqs[0:2] # first 2 sequences
      encoded_ragged_array(['LMSYAEVYGH',
                          'WKGVGKQNCAWSVNVH'], AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY*'))



    >>> my_seqs[-4:] # last 4 sequences
    encoded_ragged_array(['TGIVPMRM*S',
                          'CENVC',
                          'RSTWF',
                          'NTIFMC'], AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY*'))

Some basic properties of `EncodedRaggedArray` objects
------------------------------------------------------

    >>> my_seqs.shape # number of sequences and length of each sequence
    (10, array([10, 16, 17, 13, 14, 16, 10,  5,  5,  6]))

    >>> my_seqs.lengths # lengths of each sequence
    array([10, 16, 17, 13, 14, 16, 10,  5,  5,  6])

    >>> my_seqs.size # total number of elements (amino acid residues across all sequences in the encoded ragged array)
    112

    >>> my_seqs.encoding # the encoding used for the sequences
    AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY*')

Concatenation of `EncodedRaggedArray` objects
------------------------------------------------

    >>> np.concatenate([my_seqs, my_seqs[-2:]]).shape # concatenate two encoded ragged arrays and get the shape
    (12, array([10, 16, 17, 13, 14, 16, 10,  5,  5,  6,  5,  6]))

Getting unique elements and counting occurrences
-------------------------------------------------

    >>> bnp.count_encoded(my_seqs.get_column_values(0)) # count the number of occurrences of each amino acid at the first position (similar to numpy.unique)
	EncodedCounts(alphabet=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'], counts=array([0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 1, 2, 0, 1, 0, 0]), row_names=None)

Counting the number of occurrences of a specific element in each sequence
--------------------------------------------------------------------------

    >>> np.sum(my_seqs == "F", axis=-1) # count the number of occurrences of the amino acid "F" in each sequence
     array([0, 0, 1, 0, 0, 1, 0, 0, 1, 1])

Filtering `EncodedRaggedArray` objects based on a mask
------------------------------------------------------

    >>> mask = my_seqs.lengths < 8
    >>> short_seqs = my_seqs[mask]
    >>> short_seqs
    encoded_ragged_array(['CENVC',
                          'RSTWF',
                          'NTIFMC'], AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY*'))


Broadcasting and one-hot encoding
----------------------------------

    >>> short_seqs[1][..., np.newaxis] == "ACDEFGHIKLMNPQRSTVWY" # one-hot encoding of the second sequence
     array([[False, False, False, False, False, False, False, False, False,
            False, False, False, False, False,  True, False, False, False,
            False, False],
           [False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False,  True, False, False,
            False, False],
           [False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False,  True, False,
            False, False],
           [False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False,
             True, False],
           [False, False, False, False,  True, False, False, False, False,
            False, False, False, False, False, False, False, False, False,
            False, False]])

Reading sequences from file
===========================
Usually we get sequences directly from file. BioNumPy supports a range of file formats containing sequence data including fasta, fastq, indexed fasta and bam files.


Reading sequence entries (.fa, .fq, .bam)
-----------------------------------------
To read in a set of sequence entries, we can just use the `bnp.open` method (here we read a fastq file, but this works the same for fasta and bam files):

    >>> entries = bnp.open("example_data/reads.fq").read()
    >>> entries
    SequenceEntryWithQuality with 2 entries
                         name                 sequence                  quality
                 headerishere                 CTTGTTGA        [2 2 2 2 2 2 2 2]
                anotherheader                      CGG               [93 93 93]

We see we get all the entries in the file, with the corresponding fields. The `sequence` field here is an `EncodedRaggedArray` and thus supports numpy-like indexing etc:

    >>> (entries.sequence=="T").sum(axis=-1) # Count the number of T's
    array([4, 0])


Reading indexed files (.fa.fai)
-------------------------------
When reading a reference genome, we often can't read in the whole file (using `.read()`) and it doesn't make sense to read in the chromosome-sequences as entries in chunks (using `.read_chunks()`) but we rather want to read specific parts of the genome. In these cases we can use an index to be able to read specific regions. We then use `bnp.open_indexed` function:

    >>> reference_sequence = bnp.open_indexed("example_data/small_genome.fa")
    >>> reference_sequence["2"][10:20]
    encoded_array('atattagcca')

Functions workin on sequences
=============================

A set of functions working on sequences are gathered in the `bnp.seqeunce` module. 
