.. _introduction:

Getting started with BioNumPy in 10 minutes
--------------------------------------

By following this quick guide, you should be able to get started using BioNumPy to analyse your own biological data.

Prerequisites
==============

To follow this guide, you will need to know some Python and should also ideally know the basics of NumPy or array programming. If you have experience from R, you may already be used to concepts of array programming, such as indexing, slicing and vectorized operations, and the concepts here should be possible to understand. If you need to freshen up on NumPy/array programming, check out the `NumPy introduction <https://numpy.org/doc/stable/user/quickstart.html>`_.


Install BioNumPy
=================
If you have Python 3 already installed, BioNumPy should be straight-forward to install using pip:

.. code-block:: bash

    pip install bionumpy
    # or if you don't have the pip-command available, try:
    python3 -m pip install bionumpy

Note that when you install BioNumPy, NumPy is also automatically installed.


Read and analyse some data
=============================
In this very short tutorial, we will read a biological dataset and play around with it. We will download a fastq file from the BioNumPy example data, but feel free to use any dataset you want to play around with (see :ref:`supported file formats<supported_file_formats>`).

We start by downloading a fastq file:

.. code-block:: bash

    wget https://github.com/bionumpy/bionumpy/raw/main/example_data/big.fq.gz

In the following, we assume you have the big.fq.gz in a directory called `example_data`.

All BioNumPy programs typically start by importing both NumPy and BioNumPy. As a convention, we import BioNumPy as `bnp`:

.. testcode::

    import numpy as np
    import bionumpy as bnp

    # open the file
    f = bnp.open("example_data/big.fq.gz")
    data = f.read()  # reads the whole file into memory
    print(data)

The output from this tells us that we have a `SequenceEntryWithQuality` with 1000 entries (sequences). This is is the data type that BioNumPy automatically chooses to use when we read a fastq-file. If you instead read a vcf, bam or another format, you will get a different entry type:

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

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



The `data` variable now is a type of `dataclass` with the fields `name`, `sequence` and `quality`. Each of these behaves very much like two-dimensional NumPy-arrays (except that not all rows have the same size, since sequences can vary in length). If you for instance want to find the average base qualities, that is as simple as:

    >>> np.mean(data.quality)
    11.677166150424176

As with NumPy, you can also take the mean and other operations over various axis. Specifying `axis=0` gives you the mean over the first axis, i.e. the mean of the base quality at each base position:

    >>> np.mean(data.quality, axis=0) # doctest:+ELLIPSIS
    array([5.194, 4.599, 5.591, ..., 5.   , 6.   , 6.   ]...


BioNumPy data can also be indexed exactly as you would index NumPy arrays. This means that if you e.g. want to get all the sequence entries with more than 30% G's, you could use NumPy-syntax like this:

    >>> mask = np.mean(data.sequence == "G", axis=-1) >= 0.3
    >>> data[mask]
    SequenceEntryWithQuality with 13 entries
                         name                 sequence                  quality
      28c83abf-8f04-4651-a...  GAGCGCTGGTTCGGTTATCA...  [ 5  3  4  5  3  3  4
      666246b0-63b1-46b1-8...  CGGTGTAGCGTTTGATCTAG...  [13  3  4  7  2  5  4 1
      c1bd65c1-e3bb-40e3-8...  CGGTATGCGCTGCGTTCAGT...  [10  3  5  3  4  3  3
      8153b049-b41a-4413-a...  TAATTGCTGGATATTCCTCG...  [ 3  3  3  8 10  3  3
      ae639e55-f513-4be3-b...  CGTGTTGCGCCCGTTCAGTT...  [ 3  3  4  5  3  5  4
      30294e76-a860-4690-9...  CATTTGTACTTCCGTTCAAT...  [ 8  5  7  9  8  4  5
      5f404562-4c04-4b6d-a...  CGGTGATGCTTTGGTTACGG...  [12  3  7  8  4  2  2
      c789cd2e-01ef-4aac-8...  CTGGTGGCCGCTGGTTCGAT...  [ 7  3  2  5  6  2  3
      3d0b3924-3afc-4f48-9...  CAGTGTACTTCGTTCAGTTT...  [ 7  2  3  7  4  8  8 1
      d053ad32-9857-4440-b...  GTTGTAGCGCTACGTTTGGT...  [ 4  3  4  6  5  7  7



Final notes
============

The above examples shows how to use BioNumPy for a specific file format and dataset, but the concept is the same for all datasets. BioNumPy should be thought of as more of a toolkit rather than a collections of functions/modules. Once you learn how to use the toolkit, you should be able to very efficiently analyse many types of biological datasets.

Continue to see an overview of :ref:`what you can do with bionumpy<what_can_you_do>`.

