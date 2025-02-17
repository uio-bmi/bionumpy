

Example (for tutorial-developers)
---------




.. testsetup::

    # code that will run but not be shown
    import bionumpy as bnp
    import numpy as np
    np.random.seed(1)

Some text.

.. testcode::

    a = np.random.randint(3, 6)
    print(a)


Some text. Doctest will assert that the output from previous testcode is equal to testoutput below:

.. testoutput::

    4


Show code by using literalinclude.
This code is tested by pytest and needs to be complete (no autoimports or testsetup):

.. literalinclude:: /../scripts/example_example.py