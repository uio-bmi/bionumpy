from typing import List, Union, Iterable, Tuple, Dict
from .genomic_track import GenomicTrack, GenomicTrackNode
from ..datatypes import Interval
from .intervals import get_pileup, merge_intervals, extend_to_size, clip
from ..streams import groupby, NpDataclassStream
from ..computation_graph import StreamNode, Node, ComputationNode
from .geometry import Geometry
import dataclasses


class GenomicIntervals:
    @property
    def start(self):
        return self._intervals.start

    @property
    def stop(self):
        return self._intervals.stop

    @property
    def chromosome(self):
        return self._intervals.chromosome

    def extended_to_size(self, size: int) -> 'GenomicIntervals':
        return NotImplemented

    def merged(self, distance: int) -> 'GenomicIntervals':
        return NotImplemented

    def get_mask(self) -> GenomicTrack:
        return NotImplemented

    def get_pileup(self) -> GenomicTrack:
        return NotImplemented

    @classmethod
    def from_track(cls, track: GenomicTrack):
        return cls(track.get_data(), track._genome_context)

    @classmethod
    def from_intervals(cls, intervals: Interval, chrom_sizes: Dict[str, int]):
        return GenomicIntervalsFull(intervals, chrom_sizes)

    @classmethod
    def from_interval_stream(cls, interval_stream: Iterable[Interval], chrom_sizes: Dict[str, int]):
        interval_stream = groupby(interval_stream, 'chromosome')
        interval_stream = StreamNode(pair[1] for pair in interval_stream)
        return GenomicIntervalsStreamed(interval_stream, chrom_sizes)


class GenomicIntervalsFull(GenomicIntervals):
    def __init__(self, intervals: Interval, chrom_sizes: Dict[str, int]):
        self._intervals = intervals
        self._geometry = Geometry(chrom_sizes)
        self._chrom_sizes = chrom_sizes

    def get_data(self):
        return self._intervals

    def __len__(self):
        return len(self._intervals)

    def sorted(self):
        return NotImplemented

    def __getitem__(self, idx):
        return self.__class__(self._intervals[idx], self._chrom_sizes)

    @property
    def start(self):
        return self._intervals.start

    @property
    def stop(self):
        return self._intervals.stop

    @property
    def chromosome(self):
        return self._intervals.chromosome

    def extended_to_size(self, size: int) -> GenomicIntervals:
        return self.from_intervals(self._geometry.extend_to_size(self._intervals, size), self._chrom_sizes)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        return self.from_intervals(
            self._geometry.merge_intervals(self._intervals, distance), self._chrom_sizes)

    def get_pileup(self) -> GenomicTrack:
        """Create a GenomicTrack of how many intervals covers each position in the genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicTrack
        """
        return self._geometry.get_pileup(self._intervals)

    def clip(self) -> 'GenomicIntervalsFull':
        return self.__class__.from_intervals(self._geometry.clip(self._intervals), self._chrom_sizes)

    def __replace__(self, **kwargs):
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._chrom_sizes)


class GenomicIntervalsStreamed:
    def _get_chrom_size(self, intervals: Interval):
        return self._chrom_sizes[intervals.chromosome]

    def __str__(self):
        return 'GIS:' + str(self._intervals_node)

    def __init__(self, intervals_node: Node, chrom_sizes: Dict[str, int]):
        self._chrom_sizes = chrom_sizes
        self._start = ComputationNode(getattr, [intervals_node, 'start'])
        self._stop = ComputationNode(getattr, [intervals_node, 'stop'])
        self._strand = ComputationNode(getattr, [intervals_node, 'stop'])
        self._chromosome = ComputationNode(getattr, [intervals_node, 'chromosome'])
        self._chrom_size_node = StreamNode(iter(self._chrom_sizes.values()))
        self._intervals_node = intervals_node

    def sorted(self):
        return NotImplemented

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def chromosome(self):
        return self._chromosome

    def __getitem__(self, item):
        return self.__class__(ComputationNode(lambda x, i: x[i], [self._intervals_node, item]), self._chrom_sizes)

    def extended_to_size(self, size: int) -> GenomicIntervals:
        return self.__class__(
            ComputationNode(extend_to_size, [self._intervals_node, size, self._chrom_size_node]),
            self._chrom_sizes)
        return self.from_intervals(self._geometry.extend_to_size(self._intervals, size), self._chrom_sizes)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        return self.__class__(ComputationNode(merge_intervals, [self._intervals_node]), self._chrom_sizes)

    def get_pileup(self) -> GenomicTrack:
        """Create a GenomicTrack of how many intervals covers each position in the genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicTrack
        """
        return GenomicTrackNode(ComputationNode(get_pileup, [self._intervals_node, self._chrom_size_node]),
                                self._chrom_sizes)

        return self._geometry.get_pileup(self._intervals)

    def clip(self) -> 'GenomicIntervalsFull':
        return self.__class__(ComputationNode(clip, [self._intervals_node, self._chrom_size_node]), self._chrom_sizes)
        return self.__class__.from_intervals(self._geometry.clip(self._intervals), self._chrom_sizes)

    def __replace__(self, **kwargs):
        return self.__class__(ComputationNode(dataclasses.replace, [self._intervals_node], kwargs), self._chrom_sizes)
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._chrom_sizes)
