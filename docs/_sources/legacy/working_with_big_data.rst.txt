.. _working_with_big_data:

Working with big data
---------------------
Before following this guide we assume you have read how to read a big file under Method3 in :ref:`reading_files`.

The recommended way of working with big datasets (bigger than what you can fit in memory) in BioNumPy is to use the `read_chunks()` method when reading your data.

    >>> import bionumpy as bnp
    >>> data = bnp.open("filename").read_chunks() # doctest: +SKIP

`data` is now an iterable, and when iterating over data, BioNumPy will for every iteration read a new chunk from the file, meaning that at most one chunk is kept in memory. It is possible to specify the chunk size (in bytes) to read_chunks, e.g. `read_chunks(min_chunk_size=10000000)`.

One way of working with chunks is thus to create a for-loop iterating over the chunks. However, BioNumPy makes it possible to write code as if you only had a single chunk, and can handle the tecnhicalities of making your code work on multiple chunks.

Thus, to avoid writing for-loops and having to think about multiple chunks, BioNumPy includes utility function for many common operations, which can lead to more readable code. Examples of such operations includes taking mean, making a histogram or getting a bincount of your whole dataset. For instance, if you want to take the mean of all the base qualities across all chunks, you can simply write:

    >>> chunks = bnp.open("example_data/big.fq.gz").read_chunks()
    >>> mean = bnp.mean(chunks.quality)

In addition, BioNumPy provides a streamable decorator that lets you create a function that does something on one chunk and run that function on several chunks.

    >>> # adding @streamable lets you run this function
    >>> # on multiple chunks
    >>> @bnp.streamable()
    ... def procss_chunk(chunk):
    ...    # compute something on a single chunk
    ...    pass
    ...
    >>> results = process_chunk(chunks) # doctest: +SKIP

To test if your code works, it can be a good idea to read one chunk from your file and create a function for doing what you want on that chunk. When your code works, you can simply add the @streamable decorator and run the function on all chunks. You will then likely want to summarize the results from all the chunks in some way. Here is an example of counting the number of reads on **one chunk** (the first) in a fasta file:

    >>> def count_reads(chunk):
    ...     return len(chunk.sequence)
    ...
    >>> count_reads(bnp.open("example_data/reads.fq").read_chunk())
    2

The above code will read one chunk and count the reads. Now let's add the streamable decorator and count the reads in all chunks:

    >>> @bnp.streamable()
    ... def count_reads(chunk):
    ...     return len(chunk.sequence)
    ...
    >>> all_lengths = count_reads(bnp.open("example_data/reads.fq").read_chunks())
    >>> sum_of_lengths = sum(all_lengths)

