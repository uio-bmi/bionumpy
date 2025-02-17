.. _design_principles:

Design principles and code overview
--------------------------------------

Overview
=========

BioNumpy has a quite clear design principle which can aid you in contributing code. The functionality provided by BioNumpy is divided into three parts:

* Vectorized parsing of input files

* Data handling using bnpdataclass

* Utility functions for various kinds of data

* Encoding of data into readable


Vectorized File I/O
===================

BioNumpy already supports a wide range of file formats. The file I/O itself is handled by `NpDataClassReader` which delegates the responsibility of parsing specific file formats to `FileBuffer`. Thus supporting a new file format should be done by subclassing `FileBuffer` in one way or another. This is also the case when interpreting an existing format in a different way. For instance a .bam can be read as alignments using the normal `BamBuffer` or genomic intervals using `BamIntervalBuffer`. A file buffer has three functions it needs to support, either explicitly or by inheritance.

* `read_header` is a classmethod responsible for reading in the information contained in the header of the file. This header information is later passed to the `from_raw_buffer`, and can so inform how buffers should be parsed. At the least, this method should read through the header portion of a file such that the file pointer get set to the start of the data

* `from_raw_buffer` is a classmethod responsible for detecting where in a raw buffer the last complete entry ends, and create a `FileBuffer` from the complete entries. This method should do as little as possible, since it cannot be parallelized.

* `get_data` is responsible for converting the text in the buffer into the correct data types. The output should be a `@npdataclass` object.

* `from_data` is optional and is used for file output. It takes in a `@npdataclass` object, and should convert it to a buffer that can be written directly to file. 

A good `FileBuffer` subclass should do all time-consuming operations using vectroized numpy code. That means moving the relevant bytes/characters into arrays using numpy/RaggedArray indexing, and all encoding is done using numpy. Parsing of the header should however usually be done using as readable python code as possible.

`@bnpdataclass`
===============

The main flow of a bionumpy program is to read data from file into a `@bnpdataclass` object. `@bnpdataclass` is based on the `dataclasses.dataclass` principle where the only thing you specify to create a `bnpdataclass` class is the names and types of it's fields. The added functionality of `bnpdataclass` is that it only accepts variables that are either numpy arrays, or objects that can be indexed like numpy arrays. This way a `bnpdataclass` object can itself be indexed like a numpy array. In addition it does implicit format conversion to match what is specified in the field types.

Crucial to the concept of `bnpdataclass` is that it is a dataclass, not an object. I.e. all of it's memebers are public, but it has not methods. Thus adding functionality in a `bnpdataclass` is not allowed, and should be instead implemented as utility functions.

Utility Functions
=================

Representing sequences as encoded ragged arrays and storing them together with other data in dataclasses allows for writing a lot of functionality using standard numpy expressions, which is one goal of BioNumPy. A lot of times this can also make code that make a lot of sense in the biological sequence domain. For instance `alignments[alignments.cigar_op[:, 1] != "S"]` is written using only numpy syntax, but it is clear for someone with domain knowledge that this selects the alignments that does not start with a softclip. It will however be cases where the numpy syntax itself does not reveal the code's intent clearly to one with domain knowledge, or where the numpy code itself gets so complicated that it is worth to wrap it in a function. In these cases, utility functions should be written. For instance `(variant.ref_seq.shape.length == 1) & (variant.alt_seq.shape.length == 1)` . 

Encodings
=========

A key part of BioNumpy is representing things that are not numbers as numbers, throgugh enocdings. For instance DNA bases are represented as `[0, 1, 2, 3]` . A very important thing is that sunch encoded data should never live without teh functions to decode it. This is accomplished through the `EncodedArray` class which is a `numpy` array with an encoding attribute.
