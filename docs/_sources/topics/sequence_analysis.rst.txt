
Sequence analysis
====================

BioNumPy lets you very easily read sequence data (e.g. from a FASTA or FASTQ file) into NumPy-like data-structures:

.. testcode::

    import bionumpy as bnp
    f = bnp.open("example_data/big.fq.gz")
    chunk = f.read_chunk()

`chunk` is now an object that contains the names (header in the FASTQ file), qualities and sequences as `EncodedRaggedArray`.

.. testcode::

    print(chunk)

.. testoutput::

    SequenceEntryWithQuality with 1000 entries
                         name                 sequence                  quality
      2fa9ee19-5c51-4281-a...  CGGTAGCCAGCTGCGTTCAG...  [10  5  5 12  5  4  3
      1f9ca490-2f25-484a-8...  GATGCATACTTCGTTCGATT...  [ 5  4  5  4  6  6  5
      06936a64-6c08-40e9-8...  GTTTTGTCGCTGCGTTCAGT...  [ 3  5  6  7  7  5  4
      d6a555a1-d8dd-4e55-9...  CGTATGCTTTGAGATTCATT...  [ 2  3  4  4  4  4  6
      91ca9c6c-12fe-4255-8...  CGGTGTACTTCGTTCCAGCT...  [ 4  3  5  6  3  5  6
      4dbe5037-abe2-4176-8...  GCAGGTGATGCTTTGGTTCA...  [ 2  3  4  6  7  7  6
      df3de4e9-48ca-45fc-8...  CATGCTTCGTTGGTTACCTC...  [ 5  5  5  4  7  7  7
      bfde9b59-2f6d-48e8-8...  CTGTTGTGCGCTTCGTTCAT...  [ 8  8 10  7  8  6  3
      dbcfd59a-7a96-46a2-9...  CGATTATTTGGTTCGTTCAT...  [ 5  4  2  3  5  2  2
      a0f83c4e-4c20-4c15-b...  GTTGTACTTTACGTTTCAAT...  [ 3  5 10  6  7  6  6


It is then easy to use NumPy-functionality on this data:

.. testcode::

    print("Number of Gs in the sequences:")
    print(np.sum(chunk.sequence == "G"))
    print("Number of sequences with mean quality > 30:")
    print(np.sum(np.mean(chunk.quality, axis=1) > 10))


.. testoutput::

    Number of Gs in the sequences:
    53686
    Number of sequences with mean quality > 30:
    760



**What to read next?**

    * :ref:`Tutorial: Quality analysis of FASTQ files <_fastqc_tutorial>`
    * :ref:`API documentation on the sequence module <api_sequences>`
    * :ref:`Working with Multiple Files/Data Sources <multiple_data_sources>`






