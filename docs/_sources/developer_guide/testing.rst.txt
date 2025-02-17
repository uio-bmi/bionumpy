.. _testing:

How we do testing in BioNumPy
-------------------------------

Tests for BioNumPy are automatically run when you push to the dev branch or make a pull request, but we recommend that you also run all tests locally.

.. code-block:: bash

    make test-all

The output will tell you if anything fails.

The above command will run four different types of tests:

* Unit-tests: Quick to run, will test small units of the code)
* Property-tests: Slow to run, will test properties of the code with big sets of random data)
* Testing of examples: pytest is run on all files with pattern script/*_example.py
* Doctesting: This runs doctest on all docstrings.

When doing development, you may only want to run one of these manually to save some time (e.g. only the doctests if you write documentation, or only uni-tests if you develop some code):

Run only unit-tests:

.. code-block:: bash

    pytest
    # or:
    make test

Run doctests:

.. code-block:: bash

    cd docs_source
    make doctest

Run tests in only a single example (or file):

.. code-block:: bash

    pytest example/our_example.py



Dependencies with NpStructures
=================================

Since development on BioNumPy often is done alongside development on NpStructures, we follow these rules:

* The master branch on BioNumPy should always work with the latest published release of NpStructures
* The dev-branch and other development-branches on BioNumPy should work with the dev-branch on NpStructures.

The github-actions make sure that the correct branch is being used, but when testing locally you will need to make sure you are having the correct branch of npstructures yourself.