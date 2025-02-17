.. _making_examples:


Making examples
------------------

Overview
==========
Small and simple examples and scripts that use BioNumPy can be put in the `scripts/`-directory. If you have bigger examples (e.g. using big data, having dependencies or spread across multiple scripts) it is better to create a separate repository for the example. This repository may live in the bionumpy organization on Github.


Testing of examples
====================

For automatic testing of your example, follow these rules:

1) Use the filename `[name]_example.py`. Ending your script with `_example.py` will include it in the automatic testing.

2) Put one or more functions starting with `test_` in your example file. These tests will be run by the automatic testing.


You can run your tests like this while developing:

.. code-block:: bash

    pytest scripts/your_example.py

Using data in your examples
==============================

You may use any files in the `example_data`-directory. When referencing these files, the path should be relative to the root path of BioNumpy (e.g. `example_data/reads.fq` and not `../reads.fq`).

You can put data in this directory, but only small data files. If your example needs bigger data, your should make it in a separate repository.
