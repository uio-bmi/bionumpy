.. _benchmarking_examples:

Benchmarking Examples
---------------------

Following is the code for the examples used in the benchmarking. The scripts will usually work both as a command line script,
but also includes a test showing how to use it in python.

BAM Filtering
-------------
Read a BAM file and filter out reads with MAPQ < 60:

.. code-block:: python

        with bnp.open(output[0], 'w') as f:
            for chunk in bnp.open(input[0]).read_chunks():
                f.write(chunk[chunk.mapq == 60])


Unique Intersection
-------------------
Filter one bed file on whether the intervals overlap with another bed file:

.. literalinclude:: /../scripts/unique_intersect_example.py
    :language: python

Jaccard Index
-------------
Calculate the Jaccard index between two or more bed files:

.. literalinclude:: /../scripts/jaccard_all_vs_all_example.py
    :language: python

Count Kmers
-----------
Count the number of kmers in a fasta/fastq file:

.. literalinclude:: /../scripts/kmer_counting_example.py
    :language: python

Reverse Complement
------------------
Reverse complement a fasta/fastq file:

.. literalinclude:: /../scripts/reverse_compliment_example.py
    :language: python

Sequence Length Distribution
----------------------------
Calculate the length distribution of sequences in a fasta/fastq file:

.. literalinclude:: /../scripts/sequence_length_distribution_example.py
    :language: python

Subsampling
-----------
Subsample exactly half of the seqeuences in a fasta/fastq file. This example is complicated when working on big files as it
one have to figure out how many sequences to subsample from each chunk to get exactly half of the sequences:

.. literalinclude:: /../scripts/subsample_reads_example.py
    :language: python

Translation
-----------
Translate DNA sequences in a fasta file to protein sequences:

.. literalinclude:: /../scripts/translate_example.py
    :language: python

VCF Filtering
-------------
Filter a VCF file on the allele count:

.. literalinclude:: /../scripts/vcf_allele_frequency_filtering_example.py
    :language: python