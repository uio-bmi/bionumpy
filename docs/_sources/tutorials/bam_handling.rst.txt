Working with BAM-files
------------------------
Before following this tutorial, we assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

The following example show how to handle bam files. It checks whether base pairs pairs that are softclipped in alignments differ from the base pairs in the rest of the reads. It reads in a .bam file and uses the cigar_op attribute to find softclipped reads.


.. literalinclude:: /../scripts/bam_example.py
