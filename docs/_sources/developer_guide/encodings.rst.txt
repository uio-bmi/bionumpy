
Implementing encodings
-----------------------

When implementing a new encoding, one should subclass:

* `NumericEncoding` if your encoding is numeric. An encoding is numeric if is meant to only encode numbers, e.g. BaseQualities.
* `OneToOneEncoding` if your encoding is not numeric, but it always encodes a single byte to another byte (e.g. a letter to a number).
* `Encoding` if your encoding is not one-to-one. In such cases, the encoding is typically more complex and could e.g. be meant to encode multiple elements to a single element, as in the case of `KmerEncoding`.

Note that all NumericEncodings are one-to-one (and inherit from OneToOneEncodings).

All OneToOneEncodings only needs to implement the `_encode` and `_decode` methods, and will the automatically be able to encode and decode strings, list of strings, BaseEncoded objects, etc:

The `_encode` and `_decode` methods will always be called with raw numpy-arrays. For `_encode` the array represents the byte value of the ASCII-character to encode (e.g. A will be "65").

.. testcode::

    import bionumpy as bnp

    class MyCustomEncoding(bnp.OneToOneEncoding):
        def _encode(self, data):
            return data + 1

        def _decode(self, data):
            return data - 1

        def __repr__(self):
            return "MyCustomEncoding()"

    encoding = MyCustomEncoding()
    sequence = "ACT"
    encoded = encoding.encode(sequence)
    decoded = encoding.decode(encoded)
    print(encoded)
    print(repr(encoded))
    print(repr(encoded.raw()))
    print(decoded)
    print(repr(decoded))
    print(repr(decoded.raw()))

Note that the `encode`-method can take many different types of data, such as strings and  lists of strings. It will detect the type and call the private `_encode` method with a numpy array as input. The `_decode` method is meant to produce the byte-values of the original data. The above code will give:

.. testoutput::

    ACT
    encoded_array('ACT', MyCustomEncoding())
    array([66, 68, 85], dtype=uint8)
    ACT
    encoded_array('ACT')
    array([65, 67, 84], dtype=uint8)


Some other rules/guidelines:

* Encodings that are used should always be objects, not classes (i.e. when specifying an encoding to `as_encoded_array` or `EncodedArray`).
* `EncodedArray` and `EncodedRaggedArray` should only be used explicitely when data is already encoded correctly. No sanity checking is done when using these classes. When it is not clear whether data is encoded or not, use `as_encoded_array`. Users are not meant to be using `Encoded(Ragged)Array` directly.
* All encoded data are wrapped in `EncodedArray` or `EncodedRaggedArray` EXCEPT when something is encoded with a NumericEncoding. In such cases, the encoded data will be an `np.array` or a `RaggedArray`. This is so that numeric numpy operations can be performed directly on the data.
* Unencoded data may be "Encoded" with `BaseEncoding`. When decoding data, the result will typically be an `Encoded(Ragged)Array` with encoding `BaseEncoding`.

