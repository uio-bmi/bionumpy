Pandas Interoperability
-----------------------


The `bionumpy` package is designed to work well with `pandas`. This is useful when you want to use the powerful data manipulation tools provided by `pandas` on the data you read with `bionumpy`.

It is straightforward to convert a chunk of a file or a whole file to a Pandas DataFrame:

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

Similarily, you can convert a pandas dataframe to a BnpDataclass object:

.. testcode::

    import bionumpy as bnp
    import pandas as pd
    df = pd.DataFrame({
        "name": ["read1", "read2", "read3"],
        "sequence": ["ACGT", "TGCA", "ACGT"],
    })
    print(bnp.datatypes.SequenceEntry.from_data_frame(df))

.. testoutput::

    SequenceEntry with 3 entries
                         name                 sequence
                        read1                     ACGT
                        read2                     TGCA
                        read3                     ACGT
