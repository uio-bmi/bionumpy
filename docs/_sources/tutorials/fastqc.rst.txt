.. _fastqc_tutorial

FastQC-like quality-checking of FASTQ files
------------------------------------------------


In this tutorial we will perform simple quality checking of reads from a fastq file, similarly to what the popular tool FastQC. In addition to BioNumPy, you will also need matplotlib to do some plotting.

We assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

We start by importing all we need:

    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> import matplotlib.pyplot as plt # doctest: +SKIP


We will be using the `big.fq.gz` file in the example_data folder, but feel free to use any fastq file you like.

The first step is to read a chunk of the file:

    >>> reads = bnp.open("example_data/big.fq.gz").read_chunk()

Note that this will not give us all the reads in the file, but only the first chunk. We can now inspect the reads:

===========
GC-content
===========

We want to get the GC-content (as a ratio between 0 and 1) within each sequence in the beginning of the fastq file, and plot a histogram of these numbers.

For the first chunk we read from the file, we get the sequences as a RaggedArray where each row is a sequence. Creating a boolean mask of where we have G or C is then as simple as:

    >>> mask = (chunk.sequence == "G") | (chunk.sequence == "C") # doctest: +SKIP

`mask` is now still a ragged array with 1 where we have a C or G and 0 elsewhere.

Getting the GC-content for each read can now be done by taking the mean across the reads (last axis) of this mask:

    >>> gc_contents = np.mean(mask, axis=-1) # doctest: +SKIP


We want to create a histogram of the gc-content values from all chunks. We could call get_gc_content on each chunk, add the results to a list and create a histogram from the final list, but BioNumPy also provides a utility function for creating a histogram from the results from multiple chunks:

    >>> plt.hist(gc_contents) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP


============================
Histogram of base qualities
============================
If we want to plot a histogram of all the base qualities in a subset of reads, we can use the builtin `np.bincount` function.


    >>> reads = bnp.open("example_data/big.fq.gz").read_chunk()
    >>> base_quality_bincount = np.bincount(reads.quality.ravel())
    >>> plt.plot(base_quality_bincount) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

==============================
Average base quality per base
==============================
In the GC content histogram example, we saw that we can take the mean the rows (axis=-1). If we instead want to find the average base quality for each position in the reads, we can take the mean across the columns (axis=0).

.. code-block:: python

    reads = bnp.open("example_data/big.fq.gz").read_chunk()
    scores = bnp.mean(reads.quality, axis=0)
    print(scores)
    plt.plot(scores)
    plt.show()
