from pathlib import PurePath
from typing import Union, Optional

from .gzip_reading import gzip
# import gzip
import dataclasses
from .file_buffers import FileBuffer
from .fastq_buffer import FastQBuffer
from .multiline_buffer import MultiLineFastaBuffer
from .bam import BamBuffer, BamIntervalBuffer
from .delimited_buffers import (BedBuffer, GfaSequenceBuffer,
                                GFFBuffer, ChromosomeSizeBuffer,
                                NarrowPeakBuffer, BdgBuffer, GTFBuffer)
from .buffers.sam import SAMBuffer
from .vcf_buffers import VCFBuffer
from .pairs import PairsBuffer
from .wig import WigBuffer
from .parser import NumpyFileReader, NpBufferedWriter, NumpyBamWriter
from .exceptions import FormatException
from ..streams import NpDataclassStream
from ..bnpdataclass import BNPDataClass, bnpdataclass
from .npdataclassreader import NpDataclassReader
import logging

logger = logging.getLogger(__name__)


buffer_types = {
    ".vcf": VCFBuffer,
    ".bed": BedBuffer,
    '.bdg': BdgBuffer,
    ".narrowPeak": NarrowPeakBuffer,
    ".fasta": MultiLineFastaBuffer,
    ".fa": MultiLineFastaBuffer,
    ".fna": MultiLineFastaBuffer,
    ".faa": MultiLineFastaBuffer,
    ".fastq": FastQBuffer,
    ".fq": FastQBuffer,
    ".gfa": GfaSequenceBuffer,
    ".gff": GFFBuffer,
    ".gtf": GTFBuffer,
    ".gff3": GFFBuffer,
    ".sam": SAMBuffer,
    ".bam": BamBuffer,
    ".sizes": ChromosomeSizeBuffer,
    '.wig': WigBuffer,
    '.pairs': PairsBuffer,
    '.pa5': PairsBuffer,
}


def _get_buffered_file(
    filename, suffix, mode, is_gzip=False, buffer_type=None, **kwargs
):
    open_func = gzip.open if is_gzip else open

    if buffer_type is None:
        buffer_type = _get_buffer_type(suffix)
    writer_class = NpBufferedWriter
    if suffix == ".bam":
        writer_class = NumpyBamWriter
    if mode in ("w", "write", "wb"):
        return writer_class(open_func(filename, "wb"), buffer_type)
    elif mode in ('a', 'append', 'ab'):
        return writer_class(open_func(filename, 'ab'), buffer_type)

    # kwargs2 = {key: val for key, val in kwargs.items() if key in ["has_header"]}
    file_reader = NumpyFileReader(open_func(filename, "rb"), buffer_type) # , **kwargs2)
    if is_gzip:
        file_reader.set_prepend_mode()
    lazy= kwargs.get('lazy', None)
    return NpDataclassReader(file_reader, lazy=lazy)


def _get_buffer_type(suffix):
    if suffix in buffer_types:
        return buffer_types[suffix]
    else:
        all_buffer_types = {buffer.__name__ for buffer in buffer_types.values()}
        raise RuntimeError(f"File format {suffix} does not have a default buffer type. "
                           f"Specify buffer_type argument using  "
                           f"one of {str(['bnp.'+b for b in all_buffer_types])} or change to a valid extension: {str(list(buffer_types.keys()))}")


def bnp_open(filename: str, mode: str = None, buffer_type=None, lazy=None) -> Union[NpDataclassReader, NpBufferedWriter]:
    """Open a file according to its suffix

    Open a `NpDataclassReader` file object, that can be used to read the file,
    either in chunks or completely. Files read in chunks can be used together with
    the `@bnp.streamable` decorator to call a function on all chunks in the file
    and optionally reduce the results.

    If `mode="w"` it opens a writer object. 

    Parameters
    ----------
    filename : str
        Name of the file to open
    mode : str
        Either "w" or "r"
    buffer_type : FileBuffer
        A `FileBuffer` class to specify how the data in the file should be interpreted
    lazy : bool
        If True, the data will be read lazily, i. e. only when it is accessed. This is useful
        to speed up reading of large files, but it is more memory demanding

    Returns
    -------
    NpDataclassReader
        A file reader object

    Examples
    --------
    >>> import bionumpy as bnp
    >>> all_data = bnp.open("example_data/big.fq.gz").read()
    >>> print(all_data)
    SequenceEntryWithQuality with 1000 entries
                         name                 sequence                  quality
      2fa9ee19-5c51-4281-a...  CGGTAGCCAGCTGCGTTCAG...  [10  5  5 12  5  4  3  
      1f9ca490-2f25-484a-8...  GATGCATACTTCGTTCGATT...  [ 5  4  5  4  6  6  5  
      06936a64-6c08-40e9-8...  GTTTTGTCGCTGCGTTCAGT...  [ 3  5  6  7  7  5  4  
      d6a555a1-d8dd-4e55-9...  CGTATGCTTTGAGATTCATT...  [ 2  3  4  4  4  4  6  
      91ca9c6c-12fe-4255-8...  CGGTGTACTTCGTTCCAGCT...  [ 4  3  5  6  3  5  6  
      4dbe5037-abe2-4176-8...  GCAGGTGATGCTTTGGTTCA...  [ 2  3  4  6  7  7  6  
      df3de4e9-48ca-45fc-8...  CATGCTTCGTTGGTTACCTC...  [ 5  5  5  4  7  7  7  
      bfde9b59-2f6d-48e8-8...  CTGTTGTGCGCTTCGTTCAT...  [ 8  8 10  7  8  6  3  
      dbcfd59a-7a96-46a2-9...  CGATTATTTGGTTCGTTCAT...  [ 5  4  2  3  5  2  2  
      a0f83c4e-4c20-4c15-b...  GTTGTACTTTACGTTTCAAT...  [ 3  5 10  6  7  6  6  
        
    >>> first_chunk = bnp.open("example_data/big.fq.gz").read_chunk(300000)
    >>> print(first_chunk)
    SequenceEntryWithQuality with 511 entries
                         name                 sequence                  quality
      2fa9ee19-5c51-4281-a...  CGGTAGCCAGCTGCGTTCAG...  [10  5  5 12  5  4  3  
      1f9ca490-2f25-484a-8...  GATGCATACTTCGTTCGATT...  [ 5  4  5  4  6  6  5  
      06936a64-6c08-40e9-8...  GTTTTGTCGCTGCGTTCAGT...  [ 3  5  6  7  7  5  4  
      d6a555a1-d8dd-4e55-9...  CGTATGCTTTGAGATTCATT...  [ 2  3  4  4  4  4  6  
      91ca9c6c-12fe-4255-8...  CGGTGTACTTCGTTCCAGCT...  [ 4  3  5  6  3  5  6  
      4dbe5037-abe2-4176-8...  GCAGGTGATGCTTTGGTTCA...  [ 2  3  4  6  7  7  6  
      df3de4e9-48ca-45fc-8...  CATGCTTCGTTGGTTACCTC...  [ 5  5  5  4  7  7  7  
      bfde9b59-2f6d-48e8-8...  CTGTTGTGCGCTTCGTTCAT...  [ 8  8 10  7  8  6  3  
      dbcfd59a-7a96-46a2-9...  CGATTATTTGGTTCGTTCAT...  [ 5  4  2  3  5  2  2  
      a0f83c4e-4c20-4c15-b...  GTTGTACTTTACGTTTCAAT...  [ 3  5 10  6  7  6  6  
    
    >>> all_chunks = bnp.open("example_data/big.fq.gz").read_chunks(300000)
    
    >>> for chunk in all_chunks:
    ...       print(chunk)
    ...
    SequenceEntryWithQuality with 511 entries
                         name                 sequence                  quality
      2fa9ee19-5c51-4281-a...  CGGTAGCCAGCTGCGTTCAG...  [10  5  5 12  5  4  3  
      1f9ca490-2f25-484a-8...  GATGCATACTTCGTTCGATT...  [ 5  4  5  4  6  6  5  
      06936a64-6c08-40e9-8...  GTTTTGTCGCTGCGTTCAGT...  [ 3  5  6  7  7  5  4  
      d6a555a1-d8dd-4e55-9...  CGTATGCTTTGAGATTCATT...  [ 2  3  4  4  4  4  6  
      91ca9c6c-12fe-4255-8...  CGGTGTACTTCGTTCCAGCT...  [ 4  3  5  6  3  5  6  
      4dbe5037-abe2-4176-8...  GCAGGTGATGCTTTGGTTCA...  [ 2  3  4  6  7  7  6  
      df3de4e9-48ca-45fc-8...  CATGCTTCGTTGGTTACCTC...  [ 5  5  5  4  7  7  7  
      bfde9b59-2f6d-48e8-8...  CTGTTGTGCGCTTCGTTCAT...  [ 8  8 10  7  8  6  3  
      dbcfd59a-7a96-46a2-9...  CGATTATTTGGTTCGTTCAT...  [ 5  4  2  3  5  2  2  
      a0f83c4e-4c20-4c15-b...  GTTGTACTTTACGTTTCAAT...  [ 3  5 10  6  7  6  6  
    SequenceEntryWithQuality with 489 entries
                         name                 sequence                  quality
      5f27fb90-2cb0-43d0-a...  CGTTGCTGATTCAGCATCAA...  [ 5  3  2  3  2  2  4  
      e23294d9-0079-4345-a...  CGAGCCGCTTCGTTCCGGTT...  [ 4  5  3  3  3  4  3  
      56736851-ccc9-41a6-9...  CGGTGCCTTCGTTCATTTCT...  [ 8  3  7  7  3  1  2  
      f156362d-d380-480d-8...  CTGTTGCGCCCCGGAACAGT...  [ 7 11  9  4  4  4  3  
      300f89ef-608a-463f-8...  CATACTTTGGTTCATTCTGT...  [ 3  2  4  4  4  4  5  
      755b1702-4560-4c04-a...  GGTATACTTGCCCTACGTTC...  [10  9 13  6  3  3  4  
      98de4f6b-d094-41e8-9...  GTTGTACTTCGTTCAGTTTC...  [ 4  5  6  4  7  6  6  
      00ac3f41-f735-49e5-9...  GTTGTACTTCGTTCAGCTCT...  [ 3  4  5  4  4 10 12 1
      f92d30bc-f77f-401e-9...  GTTGTACTGCTTCGTTCAGT...  [ 6  3  4  3  6  3  2  
      7e2c14c0-0662-4cc3-8...  TGATACATTACTTCGTTCGA...  [ 3  8  4  7  2  4  3  

    """

    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix in (".gz", ".bam")
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    return _get_buffered_file(filename, suffix, mode, is_gzip=is_gzip, buffer_type=buffer_type, lazy=lazy)

    
def count_entries(filename: str, buffer_type: FileBuffer = None) -> int:
    """Count the number of entries in the file

    By default it uses the file suffix to imply the file format. But
    a specific `FileBuffer` can be provided.


    Parameters
    ----------
    filename : str
        Name of the file to count the entries of
    buffer_type : FileBuffer
        A `FileBuffer` class to specify how the data in the file should be interpreted

    Returns
    -------
    int
        The number of entries in the file

    Examples
    --------
    6

    """
    reader = NumpyFileReader
    logger.info(f"Counting entries in {filename}")
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix in (".gz", ".bam")
    #if suffix == '.bam':
    #    reader = NumpyBamReader
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    open_func = gzip.open if is_gzip else open
    if buffer_type is None:
        buffer_type = _get_buffer_type(suffix)


    file_reader = reader(open_func(filename, "rb"), buffer_type)
    if is_gzip:
        file_reader.set_prepend_mode()
    chunk_counts = (chunk.count_entries() for chunk in file_reader.read_chunks(min_chunk_size=500000))
    return sum(chunk_counts)


def read(filename: str, mode: str = None, buffer_type: Optional[FileBuffer]=None) -> BNPDataClass:
    """
    Read the content of a file
    Parameters
    ----------
    filename: str
    mode: str
    buffer_type:

    Returns
    -------
    BNPDataClass

    """
    with bnp_open(filename, mode, buffer_type) as f:
        content = f.read()
    return content
