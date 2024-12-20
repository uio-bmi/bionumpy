from functools import lru_cache
from itertools import accumulate, repeat, takewhile, chain
from typing import Union, Tuple, List, Any

import numpy as np
from npstructures.raggedshape import RaggedView, RaggedView2
from npstructures import RaggedArray, ragged_slice

from ..datatypes import Bed6, BamEntry
from ..alignments.cigar import count_reference_length, split_cigar
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..encodings import BaseEncoding
from ..encodings.alphabet_encoding import BamEncoding
from .file_buffers import FileBuffer
from ..util import cached_property


class BamBufferExtractor:
    '''
    Class to handle the extraction of data from a buffer from a BAM file.
    '''

    def __init__(self, data: np.ndarray, starts: np.ndarray, ends: np.ndarray, header_data: Any,
                 is_contigous: bool = True):
        self._data = data
        self._new_lines = starts
        self._ends = ends
        self._chromosome_names = np.array([h[0] for h in header_data])
        # self._chromosome_names = as_encoded_array([h[0] for h in header_data])
        self._header_data = header_data
        self._functions = [self._get_chromosome,
                           self._get_read_name,
                           self._get_flag,
                           self._get_position,
                           self._get_mapq,
                           self._get_cigar_symbol,
                           self._get_cigar_length,
                           self._get_sequences,
                           self._get_quality]
        self._is_contigous = is_contigous

    def __len__(self):
        return len(self._new_lines)

    def _make_contigous(self):
        assert not self._is_contigous
        lens = self._ends - self._new_lines
        new_starts = np.insert(np.cumsum(lens), 0, 0)
        self._data = RaggedArray(
            self._data, RaggedView2(self._new_lines, lens)).ravel()
        self._new_lines = new_starts[:-1]
        self._ends = new_starts[1:]
        self._is_contigous = True

    @property
    def data(self) -> np.ndarray:
        if not self._is_contigous:
            self._make_contigous()
        return self._data

    def __getitem__(self, item):
        return self.__class__(self._data, self._new_lines[item], self._ends[item], self._header_data,
                              is_contigous=False)

    def _get_ints(self, offsets, n_bytes, dtype):
        tmp = self._data[(self._new_lines + offsets)[:, None] + np.arange(n_bytes)].ravel()
        ints = (tmp).view(dtype).ravel()
        assert len(ints) == len(self._new_lines), (len(ints), offsets, len(self._new_lines), n_bytes, dtype)
        return ints

    def _get_cigar_bytes(self):
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        n_cigar_bytes = n_cigar_op * 4
        return n_cigar_bytes

    def _get_quality(self):
        return ragged_slice(self._data, self._quality_start, self._quality_start + self._get_sequence_length())

    def _get_sequences(self):
        l_seq = self._get_sequence_length()
        n_seq_bytes = (l_seq + 1) // 2
        sequences = ragged_slice(self._data, self._sequence_start, self._quality_start)
        sequences = EncodedArray(
            (((sequences.ravel()[:, None]) >> (4 * np.arange(2, dtype=np.uint8)[::-1])).ravel() & np.uint8(15)),
            BamEncoding)
        new_sequences = EncodedRaggedArray(sequences, n_seq_bytes * 2)
        view = RaggedView(new_sequences._shape.starts, l_seq)
        new_sequences = new_sequences[view]
        return new_sequences

    def _get_cigar(self):
        cigars = ragged_slice(self._data, self._cigar_start, self._sequence_start)
        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.lengths // 4)
        cigar_symbol, cigar_length = split_cigar(cigars)
        return cigar_symbol, cigar_length

    def _get_cigar_symbol(self):
        return self._get_cigar()[0]

    def _get_cigar_length(self):
        return self._get_cigar()[1]

    @cached_property
    def _read_name_start(self):
        return self._new_lines + 36

    @cached_property
    def _cigar_start(self):
        return self._read_name_start + self._get_read_name_length()

    @cached_property
    def _sequence_start(self):
        return self._cigar_start + self._get_cigar_bytes()

    @cached_property
    def _quality_start(self):
        return self._sequence_start + (self._get_sequence_length() + 1) // 2

    def _get_read_name(self):
        read_names = ragged_slice(self._data, self._read_name_start, self._cigar_start - 1)
        read_names = EncodedRaggedArray(
            EncodedArray(read_names.ravel(), BaseEncoding), read_names.shape)
        return read_names

    def _get_flag(self):
        return self._get_ints(18, 2, np.uint16)

    @lru_cache(None)
    def _get_sequence_length(self):
        return self._get_ints(20, 4, np.int32)

    def _get_read_name_length(self):
        return self._data[self._new_lines + 12]

    def _get_mapq(self):
        return self._data[self._new_lines + 13]

    def _get_position(self):
        return self._get_ints(8, 4, np.int32)

    def _get_chromosome(self):
        ref_id = self._get_ints(4, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        return chromosome

    def get_field_by_number(self, i: int) -> Union[np.ndarray, EncodedArray, EncodedRaggedArray]:
        '''
        Get the data from the field with the given number.
        Parameters
        ----------
        i: int
            The field number.

        Returns
        -------
        Union[np.ndarray, EncodedArray, EncodedRaggedArray]
            The data from the field.
        '''
        return self._functions[i]()

    @property
    def size(self) -> int:
        if self._is_contigous:
            return self._data.size
        else:
            return (self._ends - self._new_lines).sum()


class BamHeader:
    """
    Class to handle the header of a BAM file.
    """

    def __init__(self, file_object):
        self._file_object = file_object
        self._header_data = []
        self.info = self.read_header()

    def read(self, n_bytes: int) -> bytes:
        bytes = self._file_object.read(n_bytes)
        self._header_data.append(bytes)
        return bytes

    def _read_zero_term(self):
        chars = []
        while True:
            chars.append(self.read(1))
            if chars[-1] == b"\x00":
                break
        return "".join(c.decode("ascii") for c in chars[:-1])

    def read_header(self) -> List[Tuple[str, int]]:
        """
        Read the header of the BAM file. Returns a list of tuples with the reference names and lengths.

        Returns
        -------
        List[Tuple[str, int]]
            The reference names and lengths.

        """
        magic = self.read(4)
        assert magic == b"BAM\1", magic
        header_length = self._read_int()
        self.read(header_length)
        n_ref = self._read_int()
        return self._handle_refs(n_ref)

    def _handle_refs(self, n_refs):
        info = []
        for _ in range(n_refs):
            ref_n = self._read_int()
            name = self._read_zero_term()
            sequence_length = self._read_int()
            info.append((name, sequence_length))
        return info

    def _read_int(self):
        return int.from_bytes(self.read(4), byteorder="little")

    def bytes(self) -> bytes:
        """
        Get the header as bytes.

        Returns
        -------
        bytes
            The header as bytes.

        """
        return b''.join(self._header_data)


class BamBuffer(FileBuffer):
    '''
    https://samtools.github.io/hts-specs/SAMv1.pdf
    '''
    dataclass = BamEntry
    supports_modified_write = False

    def __init__(self, buffer_extractor, header_data=None):
        self._buffer_extractor = buffer_extractor
        self._header_data = header_data
        self._is_validated = True

    def __getitem__(self, idx):
        return self.__class__(self._buffer_extractor[idx], self._header_data)

    def get_field_range_as_text(self, *args):
        raise Exception('Cannot write BAM file with set values')

    @property
    def size(self) -> int:
        return self._buffer_extractor.size

    @property
    def data(self) -> np.ndarray:
        return self._buffer_extractor.data

    @property
    def n_lines(self) -> int:
        return len(self._buffer_extractor)

    @classmethod
    def _read_int(self, file_object):
        return int.from_bytes(file_object.read(4), byteorder="little")

    @classmethod
    def _read_zero_term(cls, file_object):
        chars = []
        while True:
            chars.append(file_object.read(1))
            if chars[-1] == b"\x00":
                break
        return "".join(c.decode("ascii") for c in chars[:-1])

    @classmethod
    def make_header(cls, data: BamEntry) -> bytes:
        header = data.get_context("header")
        return header.bytes()

    @classmethod
    def read_header(cls, file_object) -> BamHeader:
        return BamHeader(file_object)

    @classmethod
    def _handle_refs(cls, n_refs, file_object):
        info = []
        for _ in range(n_refs):
            ref_n = cls._read_int(file_object)
            name = cls._read_zero_term(file_object)
            sequence_length = cls._read_int(file_object)
            info.append((name, sequence_length))
        return info

    @staticmethod
    def _find_starts(chunk):
        chunk = bytes(chunk)
        new_start = lambda start, _: start + int.from_bytes(chunk[start:start + 4], byteorder="little") + 4
        _starts = accumulate(repeat(0), new_start)  # chain([0], accumulate(repeat(0), new_start))
        starts = list(takewhile(lambda start: start <= len(chunk), _starts))
        return starts

    @classmethod
    def contains_complete_entry(cls, chunks):
        return True

    @classmethod
    def from_raw_buffer(cls, chunk: np.ndarray, header_data: BamHeader) -> "BamBuffer":
        chunk = np.asarray(chunk)
        starts = np.asanyarray(cls._find_starts(chunk))
        buffer_extractor = BamBufferExtractor(chunk[:starts[-1]], starts[:-1], starts[1:], header_data.info)
        return cls(buffer_extractor, header_data)

    def get_data(self)-> BamEntry:
        """
        Get the data from the buffer.

        Returns
        -------
        BamEntry
            The data from the buffer.

        """
        return BamEntry(*(self.get_field_by_number(i) for i in range(9)))

    def get_field_by_number(self, i, dtype=None)-> Union[np.ndarray, EncodedArray, EncodedRaggedArray]:
        return self._buffer_extractor.get_field_by_number(i)

    def count_entries(self) -> int:
        return len(self._buffer_extractor)


class BamIntervalBuffer(BamBuffer):
    dataclass = Bed6

    def get_field_by_number(self, i, dtype=None):
        funcs = [
            lambda: self._buffer_extractor.get_field_by_number(0),
            lambda: self._buffer_extractor.get_field_by_number(3),
            lambda: self._buffer_extractor.get_field_by_number(3) + count_reference_length(
                *(self._buffer_extractor.get_field_by_number(i) for i in (5, 6))),
            lambda: self._buffer_extractor.get_field_by_number(1),
            lambda: self._buffer_extractor.get_field_by_number(4),
            lambda: EncodedArray(
                np.where(self._buffer_extractor.get_field_by_number(2) & np.uint16(16), ord("-"), ord("+"))[:, None],
                BaseEncoding)
        ]
        return funcs[i]()

    def get_data(self):
        return self.dataclass(*(self.get_field_by_number(i) for i in range(6)))

        chromosome = self._buffer_extractor.get_field_by_number(0)
        start = self._buffer_extractor.get_field_by_number(3)
        cigar_symbol, cigar_length = (self._buffer_extractor.get_field_by_number(i) for i in (5, 6))
        read_names = self._buffer_extractor.get_field_by_number(1)
        mapq = self._buffer_extractor.get_field_by_number(4)
        strand = self._buffer_extractor.get_field_by_number(2) & np.uint16(16)
        length = count_reference_length(cigar_symbol, cigar_length)
        return Bed6(chromosome,
                    start,
                    start + length,
                    read_names,
                    mapq,
                    strand)
