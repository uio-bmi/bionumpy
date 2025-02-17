.. _supported_file_formats:

Supported file formats
-----------------------------------

This is a list of  currently supported file formats in BioNumPy. Reading files with any of these extensions with `bnp.open` will make BioNumPy automatically detect the file type and read the data into an appropriate data structure (which will be a dataclass-like object with fields).

* vcf
* bed
* fasta / fa
* bed
* fastq / fq
* gfa (limited support only)
* gff
* gtf
* gff3
* sam / bam

=======
Example
=======
We open a bed file, read one chunk and print a description of that chunk:

.. testcode::

    import bionumpy as bnp
    from bionumpy.io.delimited_buffers import DelimitedBuffer
    data = bnp.open("example_data/test.bed")
    chunk = data.read_chunk()
    print(chunk)

The above example should work wih any of the supported file formats.

This shows us that we have a a chunk of 71 intervals, and we get to see the first of these:

.. testoutput::

    Interval with 71 entries
                   chromosome                    start                     stop
                           17                  7512371                  7512447
                           17                  7512377                  7512453
                           17                  7512393                  7512469
                           17                  7512420                  7512496
                           17                  7512422                  7512473
                           17                  7512425                  7512501
                           17                  7512428                  7512504
                           17                  7512474                  7512550
                           17                  7512537                  7512613
                           17                  7512559                  7512635


Reading a new file format
------------------------------

BioNumpy works well with popular file formats in biological data.
However, if you have a custom file format that you would like to read into BioNumpy,
you can implement a new buffer class that inherits from `DelimitedBuffer` and
specify the dataclass that you would like to use to store the data.
Here is an example of how you can implement a new buffer class for a custom file format:

We define a custom dataclass (e.g. MyCustomFormat here) that corresponds to the columns in our file format.
We then define a new buffer class (e.g. MyCustomBuffer here) that inherits from `DelimitedBuffer` and specify
the dataclass (e.g. MyCustomFormat here) that we would like to use.
We can then use the `bnp.open` function to read all the files that have similar format.

.. testcode::

    from bionumpy.io.delimited_buffers import DelimitedBuffer
    from bionumpy.bnpdataclass import bnpdataclass
    import bionumpy as bnp

    @bnpdataclass
    class MyCustomFormat:
      dna: bnp.DNAEncoding
      amino_acid: bnp.AminoAcidEncoding
      v_gene: str
      j_gene: str

    class MyCustomBuffer(DelimitedBuffer):
       dataclass = MyCustomFormat

    my_sequence_data = bnp.open(filename="example_data/airr.tsv", buffer_type=MyCustomBuffer).read()
    print(my_sequence_data)

.. testoutput::

 MyCustomFormat with 100 entries
                          dna               amino_acid                   v_gene                   j_gene
      TGCGCCACCTGGGGGGACGAGCA               CATWGDEQYF                 TRBV10-2                  TRBJ2-7
      TGTGCCAGCTCACCTACGAATTC         CASSPTNSGSNYGYTF                   TRBV18                  TRBJ1-2
      TGCGGGCCCGTAATGAACACTGA              CGPVMNTEAFF                 TRBV10-2                  TRBJ1-1
      TGTGCCAGCAGTGAAGCGCGTCC         CASSEARPARMYGYTF                  TRBV6-1                  TRBJ1-2
      TGTGCCAGCAGTAGTGGGACAGG          CASSSGTGPDQPQHF                  TRBV6-3                  TRBJ1-5
      TGTGCCAGCAACCTAGCGGGGAA          CASNLAGKNTGELFF                  TRBV6-2                  TRBJ2-2
      TGTGCCAGCAGCCAACCGGGGGG         CASSQPGGSGNYGYTF                  TRBV4-2                  TRBJ1-2
      TGCGCCAGCAGCCGCGGCCTCAG           CASSRGLREETQYF                  TRBV5-1                  TRBJ2-5
      TGTGCCAGCAGCCAAGTCTCACG        CASSQVSRQDSSYEQYF                  TRBV4-2                  TRBJ2-7
      TGTGCCAGCAGGCCGGGACAGGG     CASRPGQGAPGWEDNYGYTF                   TRBV28                  TRBJ1-2

