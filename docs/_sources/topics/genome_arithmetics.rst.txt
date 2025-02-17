Genome arithmetics
======================

With BioNumPy it is easy to work with genomic intervals.

Reading a bed-file into a NumPy-like datastructure that supports indexing, slicing, etc is as simple as:

.. testcode::

    import bionumpy as bnp
    intervals = bnp.open("example_data/ctcf.bed.gz").read()
    print(intervals)

.. testoutput::

    Interval with 44722 entries
                   chromosome                    start                     stop
                         chr3                 11004979                 11005249
                         chr2                 74482500                 74482770
                         chr8                 27952349                 27952619
                         chr7                111729970                111730240
                         chr5                 43313203                 43313473
                        chr17                 64016408                 64016678
                        chr12                 94939893                 94940163
                         chr1                 30755506                 30755776
                        chr12                121352421                121352691
                        chr10                  9305297                  9305567


See the following tutorials/guides on intervals and genome arithmetics:

    * :ref:`Tutorial: Analysing read pileup inside peaks (intervals) <similarity_measures_tutorial>`
    * :ref:`More about intervals <intervals>`
    * :ref:`API documentation on the arithmetics module <arithmetics_api>`
    * :ref:`Working with Multiple Files/Data Sources <multiple_data_sources>`


