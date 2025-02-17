BioNumPy at a glance
======================================

BioNumPy is a Python library for easy and efficient representation and analysis of biological data.
Since BioNumPy builds on the interface of NumPy, people already used to NumPy or array programming should find BioNumPy very easy to get started with.

With BioNumPy, our goal is that everyone should be able to write simple, clean code that scales well to large biological datasets.

**Installing BioNumpy**

.. code-block:: bash

    pip install bionumpy

**Analyze your biosequences like numerical vectors in NumPy:**

    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> reads = bnp.open("example_data/small.fa").read()
    >>> reads
    SequenceEntry with 3 entries
                         name                 sequence
                        read1  ACACATCACAGCTACGACGA...
                        read2  AACACTTGGGGGGGGGGGGG...
                        read3  AACTGGACTAGCGACGTACT...
    >>> gc_content = np.mean((reads.sequence == "C") | (reads.sequence == "G"))
    >>> gc_content
    0.5526315789473685

..
    BioNumpy can be used to analyse a wide range of data. Follow one of the guides below:

.. _what_can_you_do:

What can you do with BioNumpy?
----------------------------------
The main philosophy behind BioNumPy is that you should be able to efficiently read biological datasets into NumPy-like data structure, and then analyse the data using NumPy-like methods, such as indexing, slicing, broadcasting and other vectorized operations (sum, mean, etc). Since NumPy arrays are used to store the data, BioNumPy has a very low memory footprint and operations are very efficient, meaning that BioNumPy is suitable for working with large datasets and can be an alternative to using libraries and tools written in more low-level languages such as C and C++.

The core components of bionumpy are highly generic, and thus not limited to any particular types of analysis. A layer of complementary convenience functionality is though included for some major types of biosequence applications, which could be a useful starting point:

.. grid:: 2

    .. grid-item-card:: :material-regular:`library_books;3em`
        :text-align: center
        :link: topics/sequence_analysis.html

        **Sequence analysis**

        Reading and analysing DNA and protein sequences

    .. grid-item-card::  :material-regular:`location_searching;3em`
        :text-align: center
        :link: topics/kmers.html

        **Kmers**

        Analysing sequence patterns such as kmers, minimzers and motifs

.. grid:: 2

    .. grid-item-card::  :material-regular:`calculate;3em`
        :text-align: center
        :link: topics/genomic_data.html

        **Genomic Data**

        Analysing genomic data on a genome (Intervals, variants, annotations, etc)


    .. grid-item-card:: :material-regular:`hub;3em`
        :text-align: center
        :link: topics/multiomics.html

        **Multiomics**

        Combining data-sets from multiple sources/domains

..
    .. grid-item-card:: :material-regular:`rocket_launch;3em`
        :text-align: center
        :link: topics/gpu.html

        **GPU-acceleration**

        Ultra-fast sequence analysis using GPU


    .. grid-item-card::  :material-regular:`construction;3em`
        :text-align: center
        :link: topics/extending_bionumpy.html

        **Build on BioNumpy**

        Combine core functionality to support your use-case

What BioNumpy is not
---------------------
* Bionumpy is not meant to be a broad catalog of specific algorithms and data structures that are useful for the biosequence domain.
* BioNumpy also do not directly interface with the many useful data repositories in the field.
For the above purposes we refer to libraries like Biopython, which already provides a broad range of specific functionalities. The BioNumpy documentation includes many examples of how to employ BioNumPy for core data representation and processing, while interacting with Biopython for specific needs.

* BioNumpy also does not provide any tailored plotting or visualisation functionality for the biosequence domain.
Instead, the flexible data operations on bionumpy objects makes it easy to compute representations that can be visualised using generic libraries like Plotly.




Read next
----------

 * :ref:`Getting started with BioNumPy in 10 minutes<introduction>`
 * :ref:`Learn how to efficiently read large data files with BioNumPy<reading_files>`
 * :ref:`Check out examples of how to use BioNumPy efficiently <benchmarking_examples>`
 * :ref:`Check out the various tutorials<tutorials_menu>`
 * :ref:`Learn how to avoid memory issues<best_practices>`
