.. _similarity_measures_tutorial:

Computing the similarity between two BED-files
-----------------------------------------------

Calculating similarity measures between bed files is easy with BioNumPy. Bed files can be read as binary masks over a
genome, and the similarity between two bed files can be computed using for instance the Jaccard or Forbes similarity measures.

In the following example, we read two bed files (transcription factor binding sites for CTCF and ZNF263) and compute the similarity:

.. literalinclude:: /../scripts/forbes_example.py



