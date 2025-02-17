.. _api_sequences

Sequences
-----------

The sequence module of BioNumPy provides various functions for analysing sequences, such as getting kmers and minizers or computing motif scores across sequences.

Example:

.. testcode::

    import bionumpy as bnp
    file = bnp.open("example_data/big.fq.gz")
    sequence = file.read().sequence
    sequence = bnp.change_encoding(sequence, bnp.DNAEncoding)
    kmers = bnp.sequence.get_kmers(sequence, 31)
    print(kmers[0:3, 0:2])  # first three sequences, first 2 kmers

.. testoutput::

    [CGGTAGCCAGCTGCGTTCAGTATGGAAGATT, GGTAGCCAGCTGCGTTCAGTATGGAAGATTT]
    [GATGCATACTTCGTTCGATTTCGTTTCAACT, ATGCATACTTCGTTCGATTTCGTTTCAACTG]
    [GTTTTGTCGCTGCGTTCAGTTTATGGGTGCG, TTTTGTCGCTGCGTTCAGTTTATGGGTGCGG]



API documentation
===================

.. currentmodule:: bionumpy.sequence

.. autofunction:: get_kmers
.. autofunction:: get_minimizers
.. autofunction:: get_motif_scores
.. autofunction:: count_encoded
.. autofunction:: match_string
.. autoclass:: PWM
   :members:
