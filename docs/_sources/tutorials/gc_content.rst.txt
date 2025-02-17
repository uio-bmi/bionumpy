.. _gc_content:
Computing GC content inside genes
----------------------------------------------------------------------------------

This example shows how a classic bioinformatics task like computing GC content inside genes
can be performed intuitively and efficiently using BioNumPy.

It combines two main workhorses of BioNumPy: Sequence arrays for DNA, and Interval arrays for genomic locations.

Sequence arrays are read in from fasta files. Interval arrays are read in from bed files.

The crux of the example is that the gene parts of a DNA sequence are extracted by a single, intuitive indexing operation, simply: `selected_seq = sequence[intervals]`.

The count of a given nulceotide `nn` inside these gene regions is then computed simply by `(selected_seq == nn).sum()`.

Finally, to compute GC content outside genes, one first computes the regions outside genes by inverting a mask of gene regions: `~get_boolean_mask(genes, len(chr1_sequence))`

A first version of the example reads in the full sequence and interval data from fasta and bed files that only contain data for a single chromosome: `bnp.open(seq_fn).read()`:

.. literalinclude:: /../scripts/gc_one_chr_example.py

A second version still reads in the full data by the `read` function, but adds in how to handle fasta and bed files containing data for multiple chromosomes by splitting the interval data based on the groupby function, and accessing corresponding entries from the read in sequence data:

.. literalinclude:: /../scripts/gc_multiple_chr_example.py

A third version shows a more scalable version of the same code, which uses the `read_chunks` function (see :ref:`reading_files`) with the fasta file to get a generator that reads in a single chromosome sequence at a time, allowing it to handle full-genome data with limited memory footprint:

.. literalinclude:: /../scripts/gc_by_chunks_example.py

The fourth version shows how to use the `Genome` functionality to conveniently calculate the GC content:

.. literalinclude:: /../scripts/genomic_gc_example.py
