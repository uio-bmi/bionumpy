.. _reading_files:

Reading files
-------------

There are three main ways of getting your data into memory in BioNumPy. Which way you choose mainly depends on how big your data is:

Method 1: Reading all your data at once
=======================================
If you are working with datasets that are small enough to fit into memory, you can read the whole file. The benefit of this approach is that you don't need to take into account that the data has been split into chunks (next methods):

    >>> import bionumpy as bnp
    >>> data = bnp.open("example_data/test.bed").read()
    >>> print(data)
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

In the above example, `data` will be an NpDataClass object with various fields, depending on the file type. For instance, if you read a fastq file, you will be able to access the sequences (`data.sequence`) and the base qualities (`data.quality`).


Method 2: Reading a single chunk from a large file
==================================================
If you have a large file, and don't want to read all of it, you can read a single chunk. This can be useful if you only want to get to know your data or perform an analysis on a small subset of your data:

    >>> file = bnp.open("example_data/reads.fq")
    >>> chunk = file.read_chunk(min_chunk_size=10000000)

Here `min_chunk_size` is the number of bytes to read. If the file is small enough, you may get the whole file by setting chunk size to a big number, but if you want to make sure you read the whole file as one single cunk you should choose method 1 instead.


Method 3: Reading whole file as a stream of chunks
==================================================
This is the preferred way you should use for all large files (typically fasta files, fastq files and bam files). The idea is to read the whole file as chunks that can be iterated over. This way, only one single chunk will be kept in memory at the time.

    >>> file = bnp.open("example_data/big.fq.gz")
    >>> for chunk in file.read_chunks(min_chunk_size=100000):
    ...    # do stuff with this chunk, e.g. take the mean
    ...    print(chunk.quality.mean())
    11.243155401311078
    11.799580504498538
    11.447879005326635
    11.753348856321052
    11.67464738973286
    12.154069194606311

Remember that when one chunk has been processed, that chunk is "lost" and the generator continues to the next chunk. A typical implementation will thus be to write a function that does all you want to do one one chunk, run that function on each chunk and summarize the results in some way.


