Filtering FASTQ reads
----------------------

Before following this tutorial, we assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

The following is an example of a small script that filters FASTQ reads. This example illustrates the use of multiple functions decorated with `@streamable()`. Each function is designed so that it initially works on one chunk, but with the streamable descorator, we can send chunks from a file and BioNumPy handles the rest for us.

This example also illustrates how to chain multiple functions.


.. literalinclude:: /../scripts/fastq_filtering_example.py
