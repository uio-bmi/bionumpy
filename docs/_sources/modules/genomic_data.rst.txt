.. _arithmetics_api:

Genomic Data
------------

The `bnp.genomic_data` module contains functions for analysing data that belongs to a reference genome, e.g .bed, .bam, .vcf, .gtf files and so on. The main entry point is the `Genome` class that creates genomic data objects, either from file, or from `BNPDataClass` objects. The `Genome` file is most commonly instaniated from a `.chrom.sizes` or `.fasta` file, which contains information about the size of each chromosome in the genome.


API documentation
===================

.. currentmodule:: bionumpy.genomic_data		   

.. autoclass:: Genome
    :members:	       

.. autoclass:: GenomicArray
    :members:	       

.. autoclass:: GenomicIntervals
    :members:	       

.. autoclass:: GenomicLocation
    :members:	       

.. autoclass:: GenomicAnnotation
    :members:	       

.. autoclass:: GenomicSequence
    :members:	       
