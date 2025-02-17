Broadcastable Functions
=======================

A key element of using BioNumPy efficiently is being able to write broadcastable functions on sequences.
This will let NumPy and BioNumPy handle the business of applying the function to large sets of sequences.
Below is a short description of how normal broadcasting works in NumPy, and how this can be utilized to make broadcastable sequence functions.

Elementwise broadcasting in NumPy
---------------------------------
One of the main functionalities of NumPy is to apply elementwise operations to large arrays of data. For functions with more that one parameters (such as adding or multiplying) this requires that the arrays are somewhat of the same shape. Of course, if the arrays have exactly the same shape, the elementwise operations is just applied to corresponding elements of the array
  
    >>> import numpy as np
    >>> a = [[1, 2, 3], [4, 5, 6]]
    >>> b = [[1, 10, 100], [2, 20, 200]]
    >>> np.add(a, b)
    array([[  2,  12, 103],
           [  6,  25, 206]])

Broadcasting comes into play when one (or both) of the arrays is missing data in one of the dimensions

Vectorwise broadcasting
-----------------------
Elementwise broadcasting is quite common to utilize. But one can also write functions that broadcast functions on vectors. An example is the dot product of vectors. The dot product maps two vectors into a scalar by summing the elementwise products.

    >>> v = np.array([2, 3, 5])
    >>> w = np.array([3, 5, 7])
    >>> print(v.dot(w))
    56
    >>> print((v*w).sum())
    56

Now the broadcasting does not work as it would with the elementwise operations. A dot product between `[3, 2, 1]` and `[1]`, will yield an error since it is not natural to think that the ones should be broadcasted to a vector. For vector functions, the last axis needs to be of the same size. For all the other axes, the broadcasting works as usual.

To illustrate the point, we can look at the matrix multiplication `@`. Here the matrices are the fundamental elements, and so a `(2, 5, 3)` array is interpreted as two 5x3 matrices as apposed to ten  size 3 vectors.

Sequence broadcasting
---------------------
This leads us to the main point of sequence broadcasting. We want to write our sequence functions in a way that it will interpret n-dimensional arrays as (n-1)-dimensional arrays of sequences, and broadcast accordingly. Lets say we have a function that counts the number of "G"'s in a sequence. A simple solution would be::

    def count_gs(sequence):
        return np.sum(sequence == "G")

However, if we give this function an array of 10 sequences, it will just count the number of 'G's in the whole set, not for each individual sequence. I.e. it does not broadcast as we want. If however we give the sum the keyword `axis=-1`, it will always sum over each sequence individually, no matter the shape of the array::

    def count_gs(sequence):
        return np.sum(sequence == "G", axis=-1)

This trick of doing reductions along the last axis will often be enough to make sequence->scalar functions broadcast as sequences (reductions are for instance `np.sum, np.max, np.argmax, np.mean` or other functions that map arrays to scalars).

Another general rule is to add axis to the very end, when adding axes to an array. This is not so common with sequences, but one example is one-hot encoding. One hot encoding maps letters in an alphabet of size `N` to `n`-dimensional binary vectors. For a single letter, this function could be written::

    def one_hot(letter, alphabet):
        return letter == alphabet

This works for a single letter, but would fail for even a single sequence. To make it work for a single sequence, we could add a dimension to the letter. This is anyways a reasonable thing to do, since one-hot ecoding does indeed add a dimension to your input::

    def one_hot(letter, alphabet):
        return letter[:, np.newaxis] == alphabet

However, this will fail for arrays of sequences, since it assumes that the `letter` array is one-dimensional. The best, most broadcastable and general, solution is to add the last axis using a `...` to denote the exisitng axes::
  
    def one_hot(letter, alphabet):
        return letter[..., np.newaxis] == alphabet

This will work for any-dimensional arrays of sequences.

Now we have seen how to make functions that removes an axis broadcastable, and functions that add an axis. Of course, functions that keep the number of axes are often trivially broadcastable, so with these two tricks, we can often make a sequence function broadcastable :)
