IO
----
The IO module provides functions and classes for reading and writing files. Central access points are the `open` function for reading files, the `open_indexed` function for reading indexed files, and the `count_entries` function for counting entries in a file.

In addition, several `FileBuffer` classes are exposed, that can be used to specify how a file should be interpreted. Giving these as `buffer_type` argument to `open` overrides any automatic format detection based on filename suffix


API documentation
===================

.. currentmodule:: bionumpy.io

.. autofunction:: bnp_open
.. autofunction:: open_indexed
.. autofunction:: count_entries
		  
.. autoclass:: BedBuffer
.. autoclass:: VCFBuffer
.. autoclass:: FastQBuffer
.. autoclass:: TwoLineFastaBuffer
.. autoclass:: IndexedFasta
   :members: 
