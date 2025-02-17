.. highlight:: shell

================================
Installation and getting started
================================


To install BioNumPy, run this command in your terminal:

.. code-block:: console

    $ pip install bionumpy

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

===========================
Your first BioNumPy program
===========================

Assuming you have a file with som biological data, the following will read a chunk from the file and print a description of the content:

.. doctest::

    >>> import bionumpy as bnp
    >>> file = bnp.open("example_data/reads.fq")
    >>> chunk = file.read_chunk()
    >>> print(chunk)
    SequenceEntryWithQuality with 2 entries
                         name                 sequence                  quality
                 headerishere                 CTTGTTGA        [2 2 2 2 2 2 2 2]
                anotherheader                      CGG               [93 93 93]

If the above works, you are now ready to read more about :ref:`reading_files`.
