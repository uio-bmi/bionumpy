.. _multiple_data_sources:

========================================
Working with Multiple Files/Data Sources
========================================

A lot of analyses requires working on multiple data sources that are annotated on the same reference sequences. Especially when working with reference genomes. A common example is analysing some variants from a vcf file, but only those variants that fall within some intervals specified in a bed file. One way to do this in BioNumPy is to use the `MultiStream` class. This is designed to take in multiple data streams or indexed data and deliver synchronized streams that return data corresponding to the same reference sequence/chromosome in each chunk. It is best shown in an example.

    >>> import bionumpy as bnp
    >>> variants = bnp.open("example_data/few_variants.vcf").read_chunks()
    >>> intervals = bnp.open("example_data/small_interval.bed").read_chunks()
    >>> reference = bnp.open_indexed("example_data/small_genome.fa")

These variants and intervals come in chunks, and stream of chunks are not really satisfied. In addition the reference genome is read as an indexed file, and so is able to give us arbitrary parts of the reference genome at demand. We can synch these three sources up using MultiStream:

    >>> multistream = bnp.MultiStream(reference.get_contig_lengths(),
    ...                               sequence=reference,
    ...                               variants=variants,
    ...                               intervals=intervals)


The attributes we specify for the multistream are now synched streams that give data corresponding to the sequences listed in `reference.get_sequence_lengths()` one at a time. This means we can give these streams to any function with the `streamable` decorator. For example:

    >>> from bionumpy.arithmetics import get_boolean_mask
    >>> @bnp.streamable(sum)
    ... def get_letter_counts_around_variants(reference_sequence, variants, intervals):
    ...     mask = get_boolean_mask(intervals, len(reference_sequence))
    ...     positions = variants.position
    ...     positions = positions[mask[positions]] # Only look at variants inside intervals
    ...     surrounding_sequences = reference_sequence[positions-10:positions+10].ravel()
    ...     return bnp.count_encoded(bnp.as_encoded_array(surrounding_sequences, bnp.DNAEncoding))
    
    
    >>> get_letter_counts_around_variants(multistream.sequence, multistream.variants, multistream.intervals)
    EncodedCounts(alphabet=['A', 'C', 'G', 'T'], counts=array([155, 127, 116, 122]), row_names=None)


Sort Order
------------

Unindexed files (.bed, .vcf) are read from start to finish. This means that the data for each chromosome comes one at a time in the order they appear in the file. Multistream expects that this order is the same as that specified in the `seqeuence_lengths` parameter, which is unfortunately not always the case. This is usually due to the fact that different programs sorts the files in different ways: `["chr1", "chr2",,, "chr10", "chr11"]` or `["chr1", "chr10", "chr11",,, "chr2"]`. If you only have one unindexed file, the easiest solution to a sort order discrepancy is to change the sort order of the `sequence_lengths` parameter. This can be done with the `sort_dict_by_key` function, with `human_key_func` or `None` as key.

    >>> sequence_lengths = {"chr1": 10, "chr11": 20, "chr2": 40, "chr10": 50}
    >>> bnp.MultiStream.sort_dict_by_key(sequence_lengths)
    {'chr1': 10, 'chr10': 50, 'chr11': 20, 'chr2': 40}
    >>> bnp.MultiStream.sort_dict_by_key(sequence_lengths, key=bnp.MultiStream.human_key_func)
    {'chr1': 10, 'chr2': 40, 'chr10': 50, 'chr11': 20}


Two Indexed Files
------------------

If you have two unindexed files, with conflicting sort order, it is not enough to change the sort order of the `sequence_lenghts`. Hopefully, one of the files is small enough that it can fit into memory, so that we can turn it into a dict with chromosomes as key. Such a dict can be passed into the MultiStream, and is oblivious to the sort-order of the origin file. For instance, a bed file with intervals is usually quite small (unless it represents mapped reads):

    >>> intervals = bnp.open("example_data/small_interval.bed").read()
    >>> interval_dict = dict(bnp.groupby(intervals, "chromosome"))
    >>> interval_dict
    {'0': Interval with 5 entries
                   chromosome                    start                     stop
                            0                       13                       18
                            0                       37                       46
                            0                       62                       83
                            0                      105                      126
                            0                      129                      130, '1': Interval with 10 entries
                   chromosome                    start                     stop
                            1                        3                       21
                            1                       41                       65
                            1                       91                      114
                            1                      131                      153
                            1                      157                      168
                            1                      174                      201
                            1                      213                      230
                            1                      240                      268
                            1                      290                      315
                            1                      319                      339, '2': Interval with 15 entries
                   chromosome                    start                     stop
                            2                        2                       16
                            2                       44                       49
                            2                       77                      101
                            2                      108                      127
                            2                      135                      154
                            2                      163                      165
                            2                      173                      177
                            2                      201                      214
                            2                      242                      268
                            2                      292                      320, '3': Interval with 20 entries
                   chromosome                    start                     stop
                            3                        7                       34
                            3                       58                       82
                            3                       95                      101
                            3                      130                      138
                            3                      150                      170
                            3                      188                      211
                            3                      234                      261
                            3                      283                      302
                            3                      325                      352
                            3                      353                      362}
    >>> multistream = bnp.MultiStream(reference.get_contig_lengths(),
    ...                               sequence=reference,
    ...                               variants=variants,
    ...                               intervals=interval_dict)

			    
