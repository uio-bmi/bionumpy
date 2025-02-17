.. _best_practices:
--------------
Best Practices for memory management
--------------

BioNumPy relies on Numpy vectorization to make the code fast. This requires us to handle relatively large amounts of data at the same time.
This means that we need to be a bit more careful with memory management than we would be in a regular Python script.
However, by following a few simple rules, we can avoid most memory issues.

-------------------
Read data in chunks
-------------------
This is a main feature of bionumpy. For large files it is almost always better to read the data in chunks.
This is done by calling the `read_chunks` method of the `NpDataClassReader` object. This method returns a generator that yields
chunks of data, which means that we only have one chunk in memory at a time. A typical program would look like this:

.. testcode::

    import bionumpy as bnp
    f = bnp.open('example_data/reads.fq')
    with bnp.open('output.fq', 'w') as out_file:
        for chunk in f.read_chunks():
            starts_with_a = chunk.sequence[:, 0] == 'A'
        out_file.write(chunk[starts_with_a])

By reading the data in chunks, we can process files that are much larger than the available memory. We can control the size
of the chunks with the `chunk_size` argument to `read_chunks`. The default is 5000000, which is a good value for most cases, but it can be worth tuning this parameter for specific use-cases.

--------------------------------
Skip Lazy Reading when necessary
--------------------------------
By default, bionumpy parses the data lazily, which means that we only parse specific fields when they are needed. This is
done to avoid parsing the entire file when we only need a small part of it. However, this means that we have the whole text buffer for a chunk
in memory. If experiencing memory issues, it can be good idea to skip the lazy parsing and parse all the fields at once. This is done by
passing `lazy=False` to `bnp.open`. This will decrease the memory usage, but it makes file io slower.

------------------
Limit the datatype
------------------
When reading files with a lot of data per entry, for instance vcf files or bam files, it can be a good idea to limit the datatype to only the features we are interrested in.
This is either by done by using the `astype` method on the chunks, or requesting specific fields. For instance, if we only need the chromosome and position of a vcf file, we can
use the `LocationEntry` datatype, which only contains these fields. This will decrease the memory usage, in many cases so much that we can read the whole vcf file  into memory:

.. testcode::

    import bionumpy as bnp
    f = bnp.open('example_data/variants.vcf')
    locations = np.concatenate([chunk.astype(bnp.LocationEntry) for chunk in f.read_chunks()])
    print(locations)

.. testoutput::

    LocationEntry with 3 entries
                   chromosome                 position
                         chr1                   883624
                         chr1                   887559
                         chr1                   887800

