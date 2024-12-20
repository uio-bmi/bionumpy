from typing import List
from numbers import Number
import numpy.typing as npt
import dataclasses
from operator import itemgetter
import numpy as np
from npstructures import RunLength2dArray, RunLengthArray
from bionumpy.encodings.string_encodings import StringEncoding
from .bedgraph import BedGraph
from ..streams.decorators import streamable
from ..streams.grouped import chromosome_map
from ..datatypes import Interval
from ..bnpdataclass import bnpdataclass
from ..string_array import StringArray
from ..util import interleave
from ..bnpdataclass import replace


class GenomicRunLengthArray(RunLengthArray):
    def to_array(self) -> np.ndarray:
        """Convert the runlength array to a normal numpy array

        Returns
        -------
        np.ndarray
        """
        if len(self) == 0:
            return np.empty_like(self._values, shape=(0,))
        values = np.asarray(self._values)
        if values.dtype == np.float64:
            values = values.view(np.uint64)
        elif values.dtype == np.float32:
            values = values.view(np.uint32)
        elif values.dtype == np.float16:
            values = values.view(np.uint16)
        array = np.zeros_like(values, shape=len(self))
        op = np.logical_xor if array.dtype == bool else np.bitwise_xor
        diffs = op(values[:-1], values[1:])
        array[self._starts[1:]] = diffs
        array[self._starts[0]] = values[0]
        op.accumulate(array, out=array)
        return array.view(self._values.dtype)

    @classmethod
    def from_intervals(cls, starts: npt.ArrayLike, ends: npt.ArrayLike, size: int, values: npt.ArrayLike = True, default_value=0) -> 'GenomicRunLengthArray':
        """Constuct a runlength array from a set of intervals and values

        Parameters
        ----------
        starts : ArrayLike
        ends : ArrayLike
        size : int
        values : ArrayLike
        default_value :

        Returns
        -------
        'RunLengthArray'
        """
        
        assert np.all(ends > starts), (ends[ends<=starts], starts[ends <= starts])
        assert np.all(starts[1:] >= ends[:-1])
        prefix = [0] if (len(starts) == 0 or starts[0] != 0) else []
        postfix = [size] if (len(ends) == 0 or ends[-1] != size) else []

        events = np.empty(len(prefix) + len(postfix) + starts.size + ends.size, dtype=int)
        if len(prefix):
            events[0] = prefix[0]
        if len(postfix):
            events[len(prefix):-1:2] = starts
            events[len(prefix)+1:-1:2] = ends
            events[-1] = postfix[0]
        else:
            events[len(prefix)::2] = starts
            events[len(prefix)+1::2] = ends

        tmp = values
        if isinstance(values, Number):
            values = np.empty(2*(events.size//2+1), dtype=np.array(values).dtype)
            values[::2] = default_value
            values[1::2] = tmp
        else:
            values = interleave(np.broadcast(default_value, values.shape), values)
            if ends[-1] != size:
                values = np.append(values, default_value)

        if (len(starts) > 0) and (starts[0] == 0):
            values = values[1:]

        values = values[:(len(events)-1)]
        return cls(events, values, do_clean=True)

    @classmethod
    def from_bedgraph(cls, bedgraph, size=None):
        if len(bedgraph) == 0:
            assert size is not None
            return cls(np.array([0, size], dtype=int),
                       np.array([0]))
        # assert bedgraph.start[0] == 0, bedgraph
        missing_idx = np.flatnonzero(bedgraph.start[1:] != bedgraph.stop[:-1])
        if len(missing_idx):
            start = np.insert(bedgraph.start, missing_idx+1, bedgraph.stop[missing_idx])
            value = np.insert(bedgraph.value, missing_idx+1, 0)
        else:
            start, value = (bedgraph.start, bedgraph.value)
        if size is not None:
            assert bedgraph.stop[-1] <= size, (bedgraph.stop[-1], size)
        if (size is None) or (size == bedgraph.stop[-1]):
            events = np.append(start, bedgraph.stop[-1])
            values = value
        else:
            events = np.append(start, [bedgraph.stop[-1], size])
            values = np.append(value, 0)
        if events[0] != 0:
            events = np.insert(events, 0, 0)
            values = np.insert(values, 0, 0)
        return cls(events, values)

    def to_bedgraph(self, sequence_name):
        return BedGraph([sequence_name]*len(self.starts), self.starts,
                        self.ends, self.values)

    @classmethod
    def from_rle(cls, rle):
        return cls(rle._events, rle._values)

    def extract_intervals(self, intervals):
        return self[intervals]


@bnpdataclass
class RawInterval:
    start: int
    stop: int


def get_pileup(intervals: Interval, chromosome_size: int) -> GenomicRunLengthArray:
    """Get the number of intervals that overlap each position of the chromosome/contig

    This uses run length encoded arrays to handle the sparse data that
    we get from intervals.

    Parameters
    ----------
    intervals : Interval,
        Intervals on the same chromosome/contig
    chromosome_size : int
        size of the chromsome/contig

    Examples
    --------
    >>> from bionumpy.datatypes import Interval
    >>> from bionumpy.arithmetics import get_boolean_mask, get_pileup
    >>> intervals = Interval(["chr1", "chr1", "chr1"], [3, 5, 10], [8, 7, 12])
    >>> pileup = get_pileup(intervals, 20)
    >>> print(pileup)
    [0 0 0 1 1 2 2 1 0 0 1 1 0 0 0 0 0 0 0 0]

    """
    if len(intervals) == 0:
        return GenomicRunLengthArray(np.array([0, chromosome_size], dtype=int), np.array([0], dtype=int))
    rla = RunLength2dArray.from_intervals(intervals.start, intervals.stop, chromosome_size)
    return GenomicRunLengthArray.from_rle(rla.sum(axis=0))


def get_boolean_mask(intervals: Interval, chromosome_size: int):
    """Get a boolean mask representing where any inteval hits

    Uses run length encoded binary arrays to represent the areas
    covered by any interval. The mask that is returned supports numpy ufuncs, 
    so that you can run logical operations on them s.a. `& | ~` and also 
    numpy indexing so you can use it to filter positions and intervals.

    Parameters
    ----------
    intervals : Interval
        Intervals on the same chromosome/contig
    chromosome_size : int
        The size of the chromosome/contig

    Examples
    --------

    >>> intervals = Interval(["chr1", "chr1", "chr1"], [3, 5, 10], [8, 7, 12])
    >>> print(intervals)
    Interval with 3 entries
                   chromosome                    start                     stop
                         chr1                        3                        8
                         chr1                        5                        7
                         chr1                       10                       12
    >>> mask = get_boolean_mask(intervals, 20)
    >>> print(mask.astype(int))
    [0 0 0 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0 0]

    Get complement of the mask:

    >>> complement = ~mask
    >>> print(complement.astype(int))
    [1 1 1 0 0 0 0 0 1 1 0 0 1 1 1 1 1 1 1 1]

    Get the intersections (`&`) and union (`|`) of the mask and another mask

    >>> other_mask = get_boolean_mask(Interval(["chr1"], [9], [15]), 20)
    >>> intersection = mask & other_mask
    >>> print(intersection.astype(int))
    [0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]
    >>> union = mask | other_mask
    >>> print(union.astype(int))
    [0 0 0 1 1 1 1 1 0 1 1 1 1 1 1 0 0 0 0 0]

    Find wether some positions overlap the mask:
    >>> print(other_mask[intervals.start])
    [False False  True]
    """
    assert np.all(intervals.stop <= chromosome_size), (np.max(intervals.stop), chromosome_size)
    if len(intervals) == 0:
        return GenomicRunLengthArray.from_intervals(np.array([], dtype=int), np.array([], dtype=int), int(chromosome_size), default_value=False)
    merged = merge_intervals(intervals[np.argsort(intervals.start)])
    m = merged.start != merged.stop
    return GenomicRunLengthArray.from_intervals(merged.start[m], merged.stop[m],
                                                size=chromosome_size, default_value=False)
# merged.start, merged.stop, 


def human_key_func(chrom_name):
    assert chrom_name.startswith("chr"), chrom_name
    parts = chrom_name[3:].split("_", maxsplit=1)
    assert len(parts) <= 2, chrom_name
    is_numeric = 1-parts[0].isdigit()
    b = parts[0] if is_numeric else int(parts[0])
    c = parts[-1] if len(parts) == 2 else ""
    return (is_numeric, b, c)


def sort_intervals(intervals: Interval, chromosome_key_function: callable = lambda x: x, sort_order: List[str] = None) -> Interval:
    """Sort intervals on "chromosome", "start", "stop"

    Parameters
    ----------
    intervals : Interval
        Unsorted intervals

    Returns
    -------
    Interval
        Sorted intervals

    """
    if hasattr(intervals.chromosome, 'encoding') and isinstance(intervals.chromosome.encoding, StringEncoding):
        args = np.lexsort((intervals.start, intervals.chromosome))
        return intervals[args]
    if sort_order is not None:
        chromosome_key_function = {name: i for i, name in enumerate(sort_order)}.__getitem__
    s = sorted((chromosome_key_function(interval.chromosome.to_string()), interval.start, interval.stop, i)
               for i, interval in enumerate(intervals))
    indices = list(map(itemgetter(-1), s))
    return intervals[indices]

def fast_sort_intervals(intervals: Interval) -> Interval:
    if hasattr(intervals.chromosome, 'encoding') and isinstance(intervals.chromosome.encoding, StringEncoding):
        args = np.lexsort((intervals.start, intervals.chromosome))
        return intervals[args]
    if isinstance(intervals.chromosome, StringArray):
        args = np.lexsort((intervals.start, intervals.chromosome.raw()))
        return intervals[args]
    assert False, 'Fast sort intervals only works with StringEncoding and StringArray'



@chromosome_map()
def merge_intervals(intervals: Interval, distance: int = 0) -> Interval:
    """Merge a set of sorted intervals

    Merges any overlapping intervals into a single interval.
    If `distance` is specified, will merge any intervals that are closer
    than `distance` apart.

    Parameters
    ----------
    intervals : Interval
        The sorted intervals to be merged
    distance : int
        Max distance for merging

    Return
    --------
    Merged Intervals

    """
    if len(intervals) == 0:
        return intervals
    assert np.all(intervals.start[:-1] <= intervals.start[1:]), "merge_intervals requires intervals sorted on start position"
    stops = np.maximum.accumulate(intervals.stop)
    if distance > 0:
        stops += distance
    valid_start_mask = intervals.start[1:] > stops[:-1]  # intervals[:-1].stop
    start_mask = np.concatenate(([True], valid_start_mask))
    stop_mask = np.concatenate((valid_start_mask, [True]))
    new_interval = intervals[start_mask]
    new_interval.stop = stops[stop_mask]
    if distance > 0:
        new_interval.stop -= distance
    assert np.all(new_interval.start[1:] > new_interval.stop[:-1]), np.sum(new_interval.start[1:] > new_interval.stop[:-1])
    return new_interval


@chromosome_map(reduction=sum)
def count_overlap(intervals_a, intervals_b):
    starts = np.concatenate([intervals_a.start, intervals_b.start])
    stops = np.concatenate([intervals_a.stop, intervals_b.stop])
    starts.sort(kind="mergesort")
    stops.sort(kind="mergesort")
    return np.sum(np.maximum(stops[:-1]-starts[1:], 0))


@streamable()
def intersect(intervals_a, intervals_b):
    all_intervals = np.concatenate([intervals_a, intervals_b])
    all_intervals = all_intervals[np.argsort(all_intervals.start, kind="mergesort")]
    stops = np.sort(all_intervals.stop, kind="mergesort")
    mask = stops[:-1] > all_intervals.start[1:]
    result = all_intervals[1:][mask]
    result.stop = stops[:-1][mask]
    return result


@streamable()
def global_intersect(intervals_b, intervals_a):
    all_intervals = np.concatenate([intervals_a, intervals_b])
    all_intervals = all_intervals[np.lexsort((all_intervals.start, all_intervals.chromosome))]
    stops = all_intervals.stop[np.lexsort((all_intervals.stop, all_intervals.chromosome))]
    mask = stops[:-1] > all_intervals.start[1:]
    result = all_intervals[1:][mask]
    result.stop = stops[:-1][mask]
    return result


def unique_intersect(intervals_a, intervals_b, genome_size):
    genome_mask = get_boolean_mask(intervals_b, genome_size)
    entry_mask = genome_mask[intervals_a].any(axis=-1)
    return intervals_a[entry_mask]


@chromosome_map()
def extend(intervals, both=None, forward=None, reverse=None, left=None, right=None):
    directed = (forward is not None) or (reverse is not None)
    undirected = (left is not None) or (right is not None)
    assert sum([both is not None, directed, undirected]) == 1
    if both is not None:
        return dataclasses.replace(intervals, 
                                   start=intervals.start-both,
                                   stop=intervals.stop+both)
    if undirected:
        starts = interval.start-left if left is not None else interval.start
        stops = interval.stop+right if right is not None else interval.stop
        dataclasses.replace(intervals, 
                            start=starts,
                            stop=stops)
    if directed:
        if forward is None:
            forward = 0
        if reverse is None:
            reverse = 0
        return dataclasses.replace(
            intervals,
            start = np.where(interval.strand=="+",
                             intervals.start-reverse,
                             intervals.start-forward),
            stop = np.where(interval.strand=="+",
                             intervals.stop+forward,
                             intervals.start+reverse)
            )


def extend_to_size(intervals: Interval, fragment_length: int, chromosome_size: np.ndarray) -> Interval:
    """Extend/shrink intervals to match the given length. Stranded

    For intervals on the + strand, keep the start coordinate and
    adjust the stop coordinate. For - intervals keep the stop
    coordinate and adjust the start

    Parameters
    ----------
    intervals : Interval
    fragment_length : int

    Returns
    -------
    Interval
    """
    is_forward = intervals.strand.ravel() == "+"
    start = np.where(is_forward,
                     intervals.start,
                     np.maximum(intervals.stop-fragment_length, 0))
    stop = np.where(is_forward,
                    np.minimum(intervals.start+fragment_length, chromosome_size),
                    intervals.stop)

    return dataclasses.replace(
        intervals,
        start=start,
        stop=stop)


def pileup(intervals):
    chroms = np.concatenate([intervals.chromosome, intervals.chromosome])

    positions = np.concatenate((intervals.start,
                                intervals.stop))
    args = np.argsort(positions, kind="mergesort")
    values = np.where(args >= len(intervals), -1, 1)
    np.cumsum(values, out=values)
    positions = positions[args]
    intervals = np.lib.stride_tricks.sliding_window_view(positions, 2)
    mask = np.flatnonzero(intervals[:, 0] == intervals[:, 1])
    intervals = np.delete(intervals, mask, axis=0)
    values = np.delete(values, mask)
    mask = np.flatnonzero(values[1:] == values[:-1])
    values = np.delete(values, mask)
    starts = np.delete(intervals[:, 0], mask+1)
    stops = np.delete(intervals[:, 1], mask)
    return BedGraph(chroms[:values.size-1],
                    starts, stops, values[:-1])


def clip(intervals: Interval, chrom_sizes) -> Interval: 
    """Clip intervals so that all intervals are contained in their corresponding chromosome

    Parameters
    ----------
    intervals : Interval

    Returns
    -------
    Interval
    """
    return replace(
        intervals,
        start=np.maximum(0, intervals.start),
        stop=np.minimum(chrom_sizes, intervals.stop))
