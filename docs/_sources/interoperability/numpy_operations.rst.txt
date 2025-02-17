Numpy Interoperability
----------------------

Using Numpy functionality
.........................

This is a collection of short examples showing how NumPy concepts like indexing, broadcasting and reductions can be used effectively in BioNumPy:

.. code-block:: python

    import bionumpy as bnp
    import numpy as np
    import matplotlib.pyplot as plt
    
    #### Boolean indexing
    intervals = bnp.open("example_data/small_interval.bed").read()
    
    # Get intervals larger than 20
    mask = (intervals.end-intervals.start) > 20
    print(intervals[mask])
    
    
    ### slice-indexing
    # Skip first and last 10 intervals and extract every 2nd
    print(intervals[10:-10:2])
    
    
    ###  List indexing
    reference_sequence = bnp.open("example_data/small_genome.fa.fai")["1"]
    
    # Get indices of all Gs in reference
    indices = np.flatnonzero(reference_sequence[:-1] == "G")
    
    # Get all letters follwing a G
    next_letters = reference_sequence[indices+1]
    
    # Count the letters following a G
    counts = bnp.count_encoded(bnp.as_encoded_sequence_array(next_letters, bnp.DNAEncoding))
    print(counts)
    
    
    ### Broadcasting
    # Get a one_hot_encoding of a sequence
    one_hot_encoding = reference_sequence[:, np.newaxis] == "ACTG"
    print(one_hot_encoding)
    
    
    ### Reductions
    # Get reads from file
    reads = bnp.open("example_data/big.fq.gz").read()
    
    # Get proportion of G's in each sequence
    g_count = (reads.sequence == "G").mean(axis=-1)
    plt.hist(g_count); plt.show()

From NumPy to BioNumPy
......................

If you have an array of already encoded sequences, you can convert it to a BioNumPy sequence array by instantiating a 'EncodedArray' object with the array and the encoding you used. Here is an example:

.. testcode::

    import bionumpy as bnp
    import numpy as np

    original_sequence = 'ACGTACGTACGT'
    lookup = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    numpy_array = np.array([lookup[base] for base in original_sequence], dtype=np.uint8)
    encoded_array = bnp.EncodedArray(numpy_array, bnp.DNAEncoding)
    print(encoded_array)

.. testoutput::

        ACGTACGTACGT

From BioNumPy to NumPy
......................

The most common way to get a NumPy array from a EncodedArray object is to either count occurances, or creating masks of certain bases. Here is an example:

.. testcode::

    import bionumpy as bnp

    encoded_array = bnp.as_encoded_array('ACGTACGAC', bnp.DNAEncoding)
    mask = encoded_array == 'A'
    counts = bnp.count_encoded(encoded_array)
    print(mask)
    print(counts.counts)

.. testoutput::

    [ True False False False  True False False  True False]
    [3 3 2 1]

If you want to access the underlying Numpy array, you can use the 'raw()' method:

.. testcode::

    import bionumpy as bnp

    encoded_array = bnp.as_encoded_array('ACGTACGAC', bnp.DNAEncoding)
    print(encoded_array.raw())

.. testoutput::

    [0 1 2 3 0 1 2 0 1]
