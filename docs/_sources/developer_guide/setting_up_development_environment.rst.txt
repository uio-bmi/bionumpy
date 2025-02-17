
.. _development_setup:

Setting up your development environment
========================================

Follow the steps below to set up an environment that lets you do development on BioNumPy, such as making examples, writing documentation or implementing features.

Step 1 (Optional): Install development version of npstructures
****************************************************************
Development on BioNumPy is often done in parallel with development on NpStructures, meaning that new functionality on BioNumPy depends on unpublished changes on NpStructures. Thus, you will usually need the development version of NpStructures, which you can get by installing the dev branch of NpStructures:

.. code-block:: bash

    git clone git@github.com:knutdrand/npstructures.git
    cd npstructures
    git checkout dev
    pip install -e .

The `-e` installs an editable version, so that if you change anything in NpStructures those changes will have an effect.

NB: You might want to install NpStructures and BioNumPy into a virtual environment. If so, make sure to make one before starting the installation.

Step 2: Install BioNumPy locally (not from PyPi)
**************************************************

In order to develop on BioNumPy you will need to clone this repository and install BioNumPy:

.. code-block:: bash

    git clone git@github.com:bionumpy/bionumpy.git

We suggest installing a local editable copy of BioNumPy so that you get BioNumPy in your Python path. It can be a good idea to create a Python virtual environment before doing this. If you installed NpStructures in the previous step, make sure you use the same virtual environment now:

.. code-block:: bash

    cd bionumpy
    pip install -e .


Step 3: Install development dependencies
******************************************
The above installation will install all required dependencies that BioNumPy uses, but for development there are some more dependencies that you will need (e.g. for testing). Install these by running:

.. code-block:: bash

    pip install -r requirements_dev.txt


Step 4: Test your installation
********************************
BioNumPy has several different sets of test (unit tests, property testing, testing of examples and doctesting). All should pass if you've setup everything correctly. You can run all tests by running the `test`-script:

.. code-block:: bash

    ./run_tests

If the above runs without any errors, you are ready to start doing development. Continue by reading the







