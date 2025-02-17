.. _writing_documentation:

Writing documentation
=========================

To write documentation, edit or create new files in the `docs_source` directory.

To build and test the documentation locally, run:

.. code-block:: bash

    make docs

This will create html files in the `docs_source/_build` directory and should open the reuslts in your browser.

Alternatively, you can run `make html` inside the docs_source directory to have some more control.

.. code-block:: bash

    make clean


Structure of documentation
------------------------------
* `source`: Put documentation of BioNumPy concepts in this directory. You will need to manually add a link from `index.rst` to new files created.
* `modules`: This is the API-documentation. Each module (e.g. io, sequences, etc) have their own file. Each such file should briefly describe the module, and include relevant functions/classes with `.. autofunction:: functionname` or `.. autoclass:: ClassName`. Module files are automatically indexed in the menu.
* `tutorials`: Put tutorials here. These will automatically be indexed.
* `developer_guide`: Put developer guides here. You will need to manually add a link from `index.rst` to new files created.


Writing code in the documentation
----------------------------------
Code written in the documentation will be tested automatically with doctest. For instance, if writing the following, doctest will check whether the output from the Python code matches the provided output:

    >>> 2 + 2
    4

You can also write code that is tested by using the testcode/testoutput directives. `Read this guide to see how it is done <https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html#directive-testcode>`_

A thid option is to include code from e.g. the scripts directory directly into the documentation:

.. code-block::

    .. literalinclude:: /../scripts/your_example.py


Check that all code you write passes the tests by running this inside the docs_source directory:

.. code-block:: bash
    make doctest

You may get a lot of warnings, but the output should end with `Build succeeded`.
