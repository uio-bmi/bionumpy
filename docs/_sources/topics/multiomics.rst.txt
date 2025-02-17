Multiomics
===========
Molecular biology data come in range of modalities - DNA, RNA and protein sequences, as well as a broad range of annotations on such sequences, including DNA variation, DNA function, DNA metyhylation, histone modifications, protein binding to DNA, gene expression and much more. Since the philosophy of bionumpy is to provide a set of fundamental, broad-purpose data representations and operations (rather than a catalog of standard analyses), bionumpy lends itself very well to the myriad ways in which one might want to combine biosequence modalities to answer varied biological or medical questions.

Example: compute GC content inside genes
----------------------------------------
The following code computes the GC content of selected parts of a DNA sequence:
    >>> chromosome = bnp.as_encoded_array("ACGTT")
    >>> genes = bnp.datatypes.Interval.from_entry_tuples([ ("chr1",0,2), ("chr1",3,5) ])
    >>> np.mean((chromosome[genes]=="C") | (chromosome[genes]=="G"))
    0.25

And since chromosome[genes] is represented as a sequence table with one row per gene region, one can easily compute GC content per such gene region by taking mean over axis one:

    >>> np.mean((chromosome[genes]=="C") | (chromosome[genes]=="G"), axis=1)
    array([0.5, 0. ])

The above example inolved two modalities - The DNA sequence of a (very short..) chromosome and a set of  intervals on this chromosome annotated as genes. With bionumpy, DNA sequence of a genome can be directly read in from a fasta file, while regions annotated as genes for the same genome is read in from e.g. a bed or GTF file. As shown above, the DNA sequence within gene regions is easily extracted by simply indexing the genome sequence bionumpy array by the interval bionumpy array of regions annotated as genes. The following tutorial shows how to do this both at small and large scale:

:ref:`Tutorial: GC content inside genes<gc_content>`

