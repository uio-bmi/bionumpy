.. _using_bionumpy_in_your_existing_project:

Using BioNumPy in your existing project
--------------------------------------

If you already have a Python project that reads biological datasets (e.g. fastq, bam or vcf files), you can often easily add BioNumPy to your project to speed up reading and parsing of such files.

All you need is to install BioNumPy, import it, and read your file with `bnp.open`:

.. testcode::

    import numpy as np
    import bionumpy as bnp

    # open your file, bnp.open automatically detects the file format
    f = bnp.open("example_data/variants.vcf")
	# a chunk is an efficient representation of a chunk of many lines
    for chunk in f.read_chunks():
		# we can iterate over the entries
        for single_entry in chunk.to_iter():
			print(single_entry)
			# and we can access things like, chromosome, position and so on
			position = single_entry.position

.. testoutput::

    VCFEntry(chromosome='chr1', position=883624, id='rs4970378', ref_seq='A', alt_seq='G', quality='.', filter='.', info='.')
    VCFEntry(chromosome='chr1', position=887559, id='rs3748595', ref_seq='A', alt_seq='C', quality='.', filter='.', info='.')
    VCFEntry(chromosome='chr1', position=887800, id='rs3828047', ref_seq='A', alt_seq='G', quality='.', filter='.', info='.')

The above example shows how you can read any file and iterate over the lines. This way, one can in most cases easily use BioNumPy as an inplace replacement for the existing code that iterates over lines from biological data files. However, iterating over lines is not very efficient, so we can do better by working directly on chunks of the file. See :ref:`reading_files` for more information.

Using BioNumPy with Pandas
============================

Another way that BioNumPy can be used directly in existing projects is if you are using Pandas.

It is straightforward to convert a chunk of a file or a whole file to a Pandas DataFrame. This time, we read a fastq file:

.. testcode::

	import bionumpy as bnp
	f = bnp.open("example_data/big.fq.gz")
	for chunk in f.read_chunks():
		df = chunk.topandas()
		print(df)

.. testoutput::

                                  name                       sequence                        quality
    0    2fa9ee19-5c51-4281-abdd-ea...  CGGTAGCCAGCTGCGTTCAGTATGGA...  [10, 5, 5, 12, 5, 4, 3, 4,...
    1    1f9ca490-2f25-484a-8972-d6...  GATGCATACTTCGTTCGATTTCGTTT...  [5, 4, 5, 4, 6, 6, 5, 6, 1...
    2    06936a64-6c08-40e9-8a10-0f...  GTTTTGTCGCTGCGTTCAGTTTATGG...  [3, 5, 6, 7, 7, 5, 4, 3, 3...
    3    d6a555a1-d8dd-4e55-936f-ad...  CGTATGCTTTGAGATTCATTCAGGAG...  [2, 3, 4, 4, 4, 4, 6, 6, 7...
    4    91ca9c6c-12fe-4255-83cc-96...  CGGTGTACTTCGTTCCAGCTAGATTT...  [4, 3, 5, 6, 3, 5, 6, 5, 5...
    ..                             ...                            ...                            ...
    995  2eef382a-21f7-4a5b-a8d8-64...  CGTTTGCGCTGGTTCATTTTATCGGT...  [2, 6, 2, 2, 4, 2, 3, 3, 4...
    996  18949e40-d30d-49f7-8a1c-c2...  GCGTACTTCGTTCAGTTTCGGAAGTG...  [2, 2, 2, 3, 3, 3, 5, 7, 7...
    997  f4aeadf5-174e-4974-aef1-8b...  CAGTAATACTTCGTTCCAGTTCTGGG...  [9, 6, 11, 10, 2, 3, 3, 3,...
    998  6b3cb23e-3f71-435b-835f-78...  CTGTTGTACTTCGATTCATTCAGGTG...  [5, 3, 3, 5, 6, 4, 5, 3, 5...
    999  d65b5418-65d5-4bf3-aac8-aa...  CGGTGACGCTGGTTTAAATCTAACGG...  [7, 3, 4, 3, 4, 2, 2, 3, 3...

    [1000 rows x 3 columns]

Continue to see an overview of :ref:`what you can do with bionumpy<what_can_you_do>`.

Using BioNumPy RaggedArray functionalities
==========================================

Alternatively, if one is working with biological sequences and are already loaded into a list of strings,
one can use the `as_encoded_array` function to convert the list of strings into a RaggedArray.
This is useful if you want to make use of all the efficient ways of BioNumpy's functionalities on
handling non-uniform length sequences and computations on them. Below is an example. See all the available encodings in :ref:`encodings`.

.. testcode::

    import bionumpy as bnp
    my_sequences = ["TGTGCCAGCAGCGGGGATCGTAATCAGCCCCAGCATTTT",
                    "TGCAGCGTCAAGGTCCAAGCTTTCTTT",
                    "TGTGCCACCAGTGATTATTATTGGTACGAGCAGTACTTC"]
    my_ragged_array = bnp.as_encoded_array(my_sequences, bnp.encodings.alphabet_encoding.DNAEncoding)
    print(my_ragged_array)
    print(my_ragged_array.shape)

.. testoutput::

    TGTGCCAGCAGCGGGGATCGTAATCAGCCCCAGCATTTT
    TGCAGCGTCAAGGTCCAAGCTTTCTTT
    TGTGCCACCAGTGATTATTATTGGTACGAGCAGTACTTC
    (3, array([39, 27, 39]))

