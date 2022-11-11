import pytest

from bionumpy.streams.multistream import MultiStream, SequenceSizes, StreamError, SequenceInfo
from bionumpy.streams import BnpStream, NpDataclassStream
from bionumpy.bnpdataclass import bnpdataclass


@bnpdataclass
class SimpleClass:
    chromosome: str


@pytest.fixture
def stream():
    return NpDataclassStream([SimpleClass(["chr1"]*3+["chr2"]*2+["chr3"])],
                             dataclass=SimpleClass)


@pytest.fixture
def indexed():
    return {"chr1": 10, "chr2": 20, "chr3": 30}


@pytest.fixture
def sequence_sizes():
    return SequenceSizes([("chr1", 20), ("chr2", 30), ("chr3", 40)])


def test_multistream(stream, indexed, sequence_sizes):
    multistream = MultiStream(sequence_sizes, names=stream, values=indexed)
    output = list(zip(multistream.infos,
                      multistream.names,
                      multistream.values))
    true = [(SequenceInfo("chr1", 20), SimpleClass(["chr1"]*3), 10),
            (SequenceInfo("chr2", 30), SimpleClass(["chr2"]*2), 20),
            (SequenceInfo("chr3", 40), SimpleClass(["chr3"]*1), 30)]
    assert len(output)==len(true)
    for o, t in zip(output, true):
        assert o == t


def test_multistream_wtih_repr(stream, indexed, sequence_sizes):
    multistream = MultiStream(sequence_sizes, names=stream, values=indexed)
    assert "SimpleClass" in repr(multistream.names)
    assert "chr1" in repr(multistream.values)
    output = list(zip(multistream.infos,
                      multistream.names,
                      multistream.values))
    true = [(SequenceInfo("chr1", 20), SimpleClass(["chr1"]*3), 10),
            (SequenceInfo("chr2", 30), SimpleClass(["chr2"]*2), 20),
            (SequenceInfo("chr3", 40), SimpleClass(["chr3"]*1), 30)]
    assert len(output)==len(true)
    for o, t in zip(output, true):
        assert o == t


def test_raises_on_missing_seq_len(stream, indexed, sequence_sizes):
    del sequence_sizes["chr2"]
    print(sequence_sizes)
    multistream = MultiStream(sequence_sizes, names=stream, values=indexed)
    with pytest.raises(StreamError):
        print(list(zip(multistream.infos,
                       multistream.names,
                       multistream.values)))


def _test_raises_on_missing_stream(stream, indexed, sequence_sizes):
    stream = NpDataclassStream((elem[:-1] for elem in stream), dataclass=SimpleClass)
    multistream = MultiStream(sequence_sizes, names=stream, values=indexed)
    with pytest.raises(StreamError):
        print(list(zip(multistream.infos,
                       multistream.names,
                       multistream.values)))


def test_raises_on_sortorder(stream, indexed, sequence_sizes):
    stream = NpDataclassStream((elem[::-1] for elem in stream), dataclass=SimpleClass)
    multistream = MultiStream(sequence_sizes, names=stream, values=indexed)
    with pytest.raises(StreamError):
        print(list(multistream.names))
