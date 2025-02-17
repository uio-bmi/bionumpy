.. _topic_kmers


Using BioNumPy to analyse kmers, minimzers and motifs
========================================================

BioNumPy has functionality for easily finding kmers and minimizers in a set of sequences, and does so very efficiently.

For instance, finding all 3-mers and counting the number of `ACT` in some sequences is as simple as:

    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are stored numerically internally, and how they are computed will depend on the encoding of the sequences.
If your sequences is encoded with an encoding with alphabet size **4** (e.g. `bnp.DNAEncoding`)
`get_kmers` is optimized to run a lot faster by performing bit-operations.


Analysing kmers and minimizers in big datasets
-------------------------------------------------
BioNumPy is able to quite efficiently find all kmers or minimizers in large sequence datasets. Since the resulting kmers are stored numerically in NumPy arrays, NumPy functionality can be used to do analysis.

Below in an example of finding all 31-mers in a fastq file:

.. testcode::

    file = bnp.open("example_data/big.fq.gz")
    # read the file chunk by chunk to keep memory low:
    for chunk in file.read_chunks():
        sequences = chunk.sequence
        # change encoding to a DNAEncoding, this works as long as the
        # sequences only contains ACGT, and makes get_kmers extremely efficient
        sequences = bnp.change_encoding(sequences, bnp.DNAEncoding)
        print("Kmers:")
        kmers = bnp.get_kmers(sequences, k=31)
        print(kmers[0:3, 0:2])

        # kmers is an EncodedRaggedArray, one row for each read
        # you can get all the raw numeric kmers efficiently like this:
        print("Raw kmers:")
        numeric_kmers = kmers.raw()
        print(numeric_kmers[0:3, 0:2])

        # and if you don't care about the RaggedStructure, you can do ravel:
        print("Flat raw kmers:")
        numeric_flat_kmers = numeric_kmers.ravel()
        print(numeric_flat_kmers[0:4])

Output from the code above:

.. testoutput::

    Kmers:
    [CGGTAGCCAGCTGCGTTCAGTATGGAAGATT, GGTAGCCAGCTGCGTTCAGTATGGAAGATTT]
    [GATGCATACTTCGTTCGATTTCGTTTCAACT, ATGCATACTTCGTTCGATTTCGTTTCAACTG]
    [GTTTTGTCGCTGCGTTCAGTTTATGGGTGCG, TTTTGTCGCTGCGTTCAGTTTATGGGTGCGG]
    Raw kmers:
    [4360244785522956521 4548825710201280058]
    [3755975642940518834 3244836919948823660]
    [2804282287455632382 3006913581077602047]
    Flat raw kmers:
    [4360244785522956521 4548825710201280058 3443049436764013966
      860762359191003491]

