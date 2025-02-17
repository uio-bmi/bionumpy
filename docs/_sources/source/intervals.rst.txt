.. _intervals:

===========
 Intervals
===========

Working with intervals is a central part of many bioinformatics analyses. In BioNumPy, an interval is any object (`bnpdataclass`) having a `start` and `stop` attribute and optionally a `strand` attribute. We will typically get intervals by reading them from files such as `.bed`, `.gtf` or `.narrowPeaks`. Or by converting alignments to reference intervals by `alignment.to_interval`. The functionality in BioNumPy for handling intervals can be divided into two operating on intervals as independent objects, or operations with some sort of dependencies:


Independent Operations
======================

Independent operations covers any thing we want to do with a set of intervals that can be done one interval at a time. This covers simple geometric operations such as shifting and resizing; filtering based on properties of intervals. Since both the start and stop (and strand) attribute are NumPy arrays, we can implement a lot of functionality simply by using NumPy indexing and ufuncs:

    >>> import bionumpy as bnp
    >>> intervals = bnp.open("example_data/small_interval.bed").read()
    >>> extended_right = bnp.replace(intervals, stop=intervals.stop+10)
    >>> shifted = bnp.replace(intervals, start=intervals.start+5, stop=intervals.stop+5)
    >>> small = intervals[(intervals.stop-intervals.start)<50]

Even though it is possible to write these kinds of functions using simple NumPy functionality, it is often better to wrap them in functions with instructive names. In the `bnp.interval` module, a set of utility funcitons ffor working with intervals is provided.

Operations with dependencies
============================
Some queries cannot be answered just by looking at intervals as independent entries. Consider the query "all regions covered by at least two intervals". In order to answer such questions, we need to look at the intervals as a coherent set. In BioNumPy, such queries are done by reducing all the relevant information from the intervals into a  `npstructures.RunLengthArray`. The two main ways of doing this is to do make a pileup, i.e. counting how many intevals overlap with each position in the genome, or to make a boolean mask indicating which positions are covered by any interval.

    >>> from bionumpy.datatypes import Interval
    >>> from bionumpy.arithmetics import get_boolean_mask, get_pileup
    >>> intervals = Interval(["chr1"]*3, [10, 15, 20], [13, 22, 21])
    >>> get_boolean_mask(intervals, 23)
    array([False, False, False, False, False, False, False, False, False,
              False,  True,  True,  True, False, False,  True,  True,  True,
               True,  True,  True,  True, False])
    >>> get_pileup(intervals, 23)
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 2,
              1, 0])

These run length encoded arrays stores the data in a smart way, but are supposed to work a much like a NumPy array as possible. So we can do NumPy-like operations on them to achieve common goals:

Working with boolean masks
--------------------------
The boolean mask representation is good to use when you want to filter a set of intervals or positions based on overlap with another set of intervals. We will examplify that here with a use of: a set of genes,  a set of aligned reads and a set of variants.

    >>> from bionumpy.datatypes import Bed6, Variant
    >>> genes = Bed6(["chr1"]*3, [10, 15, 20], [13, 22, 21], ["gene_a", "gene_b", "gene_c"], ["."]*3, ["+", "-", "+"])
    >>> genes
    Bed6 with 3 entries
                   chromosome                    start                     stop                     name                    score                   strand
                         chr1                       10                       13                   gene_a                        .                        +
                         chr1                       15                       22                   gene_b                        .                        -
                         chr1                       20                       21                   gene_c                        .                        +
    >>> reads = Bed6(["chr1"]*5, [2, 3, 5, 7, 11], [4, 6, 8,10,12], [f"read{i}" for i in range(5)], ["."]*5, ["+", "-", "+", "-", "+"])
    >>> reads
    Bed6 with 5 entries
                   chromosome                    start                     stop                     name                    score                   strand
                         chr1                        2                        4                    read0                        .                        +
                         chr1                        3                        6                    read1                        .                        -
                         chr1                        5                        8                    read2                        .                        +
                         chr1                        7                       10                    read3                        .                        -
                         chr1                       11                       12                    read4                        .                        +
    >>> variants = Variant(["chr1"]*4, [2, 4, 11, 16], ["A", "C", "G", "T"], ["G", "C", "T", "A"])
    >>> variants
    Variant with 4 entries
                   chromosome                 position                  ref_seq                  alt_seq
                         chr1                        2                        A                        G
                         chr1                        4                        C                        C
                         chr1                       11                        G                        T
                         chr1                       16                        T                        A

Let's say we now want to get all the variants that lie within a gene, and is also covered by at least one read. We would then make a boolean mask for both the genes and the reads:

    >>> gene_mask = get_boolean_mask(genes, 23)
    >>> read_mask = get_boolean_mask(reads, 23)

We can now combine these together using the `and` operator to get the mask where both are covered:

    >>> gene_and_read_mask = gene_mask & read_mask

Finally, we want to find all the variants that lie within this mask:

    >>> covered_variants = variants[gene_and_read_mask[variants.position]]
    >>> covered_variants
    Variant with 1 entries
                   chromosome                 position                  ref_seq                  alt_seq
                         chr1                       11                        G                        T

We see that only one variant is covered both by a gene and a read.


Genomic Intervals
=================
When working with genomic intervals, we often want to only deal intervals from one chromsome at the time. In order to do this, we can use the `bnp.groupby` function on interval entries.

    >>> import bionumpy as bnp
    >>> intervals = bnp.open("example_data/small_interval.bed").read()
    >>> print(intervals)
    Interval with 50 entries
                   chromosome                    start                     stop
                            0                       13                       18
                            0                       37                       46
                            0                       62                       83
                            0                      105                      126
                            0                      129                      130
                            1                        3                       21
                            1                       41                       65
                            1                       91                      114
                            1                      131                      153
                            1                      157                      168
    >>> for chromosome, data in bnp.groupby(intervals, "chromosome"):
    ...     print(f"---Chromosome: {chromosome}---")
    ...     print(data)
    ---Chromosome: 0---
    Interval with 5 entries
                   chromosome                    start                     stop
                            0                       13                       18
                            0                       37                       46
                            0                       62                       83
                            0                      105                      126
                            0                      129                      130
    ---Chromosome: 1---
    Interval with 10 entries
                   chromosome                    start                     stop
                            1                        3                       21
                            1                       41                       65
                            1                       91                      114
                            1                      131                      153
                            1                      157                      168
                            1                      174                      201
                            1                      213                      230
                            1                      240                      268
                            1                      290                      315
                            1                      319                      339
    ---Chromosome: 2---
    Interval with 15 entries
                   chromosome                    start                     stop
                            2                        2                       16
                            2                       44                       49
                            2                       77                      101
                            2                      108                      127
                            2                      135                      154
                            2                      163                      165
                            2                      173                      177
                            2                      201                      214
                            2                      242                      268
                            2                      292                      320
    ---Chromosome: 3---
    Interval with 20 entries
                   chromosome                    start                     stop
                            3                        7                       34
                            3                       58                       82
                            3                       95                      101
                            3                      130                      138
                            3                      150                      170
                            3                      188                      211
                            3                      234                      261
                            3                      283                      302
                            3                      325                      352
                            3                      353                      362
