import dataclasses
import os
import warnings
from pathlib import Path
from typing import Iterable, Union, List

import numpy as np
from npstructures import RaggedArray

from ..datatypes import BamEntry, Interval
from ..encoded_array import  EncodedRaggedArray, EncodedArray, BaseEncoding
from ..arithmetics.intervals import fast_sort_intervals
from ..encodings import CigarOpEncoding


def alignments_to_bam_entries(alignments: Iterable['pysam.AlignedSegment'], none_on_empty=False) -> Union[
    BamEntry, None]:
    cols = alignments_to_cols(alignments)
    return cols.to_entry()


@dataclasses.dataclass
class BamAccumulator:
    reference_name: List[str]
    query_name: List[str]
    flag: List[int]
    reference_start: List[int]
    mapping_quality: List[int]
    cigar_ops_data: List[int]
    cigar_ops_lengths: List[int]
    cigar_len_data: List[int]
    cigar_len_lengths: List[int]
    seq_data: List[str]
    seq_lengths: List[int]
    qual_data: List[int]
    qual_lengths: List[int]

    def add_alignment(self, read: 'pysam.AlignedSegment'):
        accumulator = self
        accumulator.reference_name.append(read.reference_name)
        accumulator.query_name.append(read.query_name)
        accumulator.flag.append(read.flag)
        accumulator.reference_start.append(read.reference_start)
        accumulator.mapping_quality.append(read.mapping_quality)
        cigar_ops, cigar_lengths = zip(*read.cigartuples)
        accumulator.cigar_ops_data.extend(cigar_ops)
        accumulator.cigar_ops_lengths.append(len(cigar_ops))
        accumulator.cigar_len_data.extend(cigar_lengths)
        accumulator.cigar_len_lengths.append(len(cigar_lengths))
        seq = read.seq
        accumulator.seq_data.append(seq)
        accumulator.seq_lengths.append(len(seq))
        qual = read.qual
        accumulator.qual_data.append(qual)
        accumulator.qual_lengths.append(len(qual))

    def to_entry(self):
        if not len(self.reference_name):
            return BamEntry.empty()
        seq_data = np.frombuffer(bytes(''.join(self.seq_data), 'utf-8'), dtype=np.uint8)
        qual_data = np.frombuffer(bytes(''.join(self.qual_data), 'utf-8'), dtype=np.uint8)
        return BamEntry(
            np.array(self.reference_name),
            np.array(self.query_name),
            np.array(self.flag),
            np.array(self.reference_start),
            np.array(self.mapping_quality),
            EncodedRaggedArray(EncodedArray(
                np.array(self.cigar_ops_data), CigarOpEncoding),
                np.array(self.cigar_ops_lengths)),
            RaggedArray(np.array(self.cigar_len_data), np.array(self.cigar_len_lengths)),
            EncodedRaggedArray(EncodedArray(seq_data, BaseEncoding),
                               np.array(self.seq_lengths)),
            RaggedArray(qual_data, np.array(self.qual_lengths)))


def alignments_to_cols(alignments, min_start=0, accumulator=None):
    accumulator = accumulator or BamAccumulator(*(list() for _ in range(13)))
    for read in alignments:
        if read.reference_start < min_start:
            continue
        accumulator.add_alignment(read)
    return accumulator


class IndexedBamFile:
    '''
    A wrapper class around pysam to extract all reads that overlap with any interval in a set of intervals.

    Examples
    --------
    >>> import bionumpy as bnp
    >>> bn = bnp.open_indexed('example_data/ctcf_chr21-22.bam')
    >>> intervals = bnp.open('example_data/ctcf.bed.gz').read()
    >>> bn[intervals]
    BamEntry with 12649 entries
                   chromosome                     name                     flag                 position                     mapq                 cigar_op             cigar_length                 sequence                  quality
                        chr21  SOLEXA-1GA-2:2:39:563:1                        0                 10403175                       37                        M                     [36]  AGGCGGAGCCCTAGGGACAGGAG  [96 97 96 97 97 96 96 9
                        chr21  SOLEXA-1GA-2:2:51:1257:                        0                 10403191                       37                        M                     [36]  ACAGGAGGAGGGGAGTTGCGCAC  [96 97 97 96 96 96 97 9
                        chr21  SOLEXA-1GA-2:2:90:233:6                       16                 13980514                       37                        M                     [36]  ACACCCTCCCCTCGCCGCTGCAG  [66 92 90 90 94 92 79 7
                        chr21  SOLEXA-1GA-2:2:62:293:1                       16                 13980528                       37                        M                     [36]  CCGCTGCAGTGTAGAAACCCAAT  [89 95 93 93 96 94 97 9
                        chr21  SOLEXA-1GA-1:1:49:718:1                        0                 13980531                       37                        M                     [36]  CTGCAGTGTAGAAACCCAATAGC  [97 97 97 98 97 97 96 9
                        chr21  SOLEXA-1GA-2:2:57:1221:                       16                 13980533                       37                        M                     [36]  GCAGTGTAGAAACCCAATAGCGT  [97 97 97 95 93 96 93 9
                        chr21  SOLEXA-1GA-1:1:57:1445:                       16                 13980536                       37                        M                     [36]  GTGTAGAAACCCAATAGCGTCCC  [96 92 94 96 93 97 93 9
                        chr21  SOLEXA-1GA-2:2:64:1358:                        0                 14120164                       37                        M                     [36]  ACCCTTAAAAGACCCAGATGTTG  [97 98 96 97 97 98 97 9
                        chr21  SOLEXA-1GA-1:1:63:383:1                        0                 14120199                       37                        M                     [36]  ATGGAAGCAGCTTCATATCCAAG  [97 97 95 97 98 97 97 9
                        chr21  SOLEXA-1GA-1:1:111:87:1                        0                 14120203                       37                        M                     [36]  AAGCAGCTTCATATCCAAGGGTG  [97 97 95 97 98 97 98 9
    '''

    def __init__(self, filename: str, create_index=False):
        try:
            import pysam
        except ImportError:
            raise ImportError('Please install pysam to use IndexedBamFile')
        warnings.warn('Indexed bam files are experimental and may not work as expected, use at your own risk')
        if create_index:
            index_filename= Path(filename).with_suffix('.bam.bai')
            if not os.path.isfile(index_filename):
                pysam.index(str(filename))
        self._samfile = pysam.AlignmentFile(filename, 'rb')

    def _fetch_from_sorted_intervals(self, intervals):
        cur_chromosome, last_stop, last_start = None, None, None
        accumulator = BamAccumulator(*(list() for _ in range(13)))
        for interval in intervals:
            if interval.chromosome == cur_chromosome:
                assert interval.start >= last_start, f'Intervals must be sorted {(interval.start, last_start)}'
            else:
                last_stop = 0
            alignments_to_cols(self._samfile.fetch(str(interval.chromosome),
                                                   start=int(interval.start),
                                                   stop=int(interval.stop)),
                               min_start=last_stop, accumulator=accumulator)

            cur_chromosome, last_stop, last_start = interval.chromosome, interval.stop, interval.start
            continue

        return accumulator.to_entry()

    def get_all_overlapping(self, intervals: Interval) -> BamEntry:
        '''
        Extract all reads that overlap with any interval in the input item.
        Parameters
        ----------
        intervals

        Returns
        -------

        '''
        item = fast_sort_intervals(intervals)
        return self._fetch_from_sorted_intervals(item)

    def __getitem__(self, item: Interval):
        '''
        Extract all reads that overlap with any interval in the input item.
        Parameters
        ----------
        item

        Returns
        -------

        '''
        return self.get_all_overlapping(item)
