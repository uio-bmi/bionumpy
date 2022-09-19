import numpy as np
from .npdataclassstream import NpDataclassStream
import logging
from npstructures import npdataclass

logger = logging.getLogger(__name__)


def wrapper(x):
    return x


def repr_bytes(n):
    if n < 10 ** 4:
        return str(n) + "b"
    elif n < 10 ** 7:
        return str(n // 1000) + "kb"
    elif n < 10 ** 11:
        return str(n // 1000000) + "Mb"
    return str(n // 1000000000) + "Gb"


class NumpyFileReader:
    """
    Class that handles the main task of reading in chunks of data
    from file that corresponds to full entries.
    """
    def __init__(self, file_obj, buffer_type, chunk_size=5000000, has_header=False):

        self._file_obj = file_obj
        """Inits the file reader with an opened file object, and a class of the buffer type used to parse it

        Parameters
        ----------
        file_obj : file
            Opened file object
        buffer_type : FileBuffer
            Buffer class used to parse the data
        chunk_size : int
            number of bytes to read
        has_header : bool
            whether or not the file has a header

        Examples
        --------
        FIXME: Add docs.

        """
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type
        self._has_header = has_header
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )

    def __enter__(self):
        return self

    def __exit__(self):
        self._file_obj.close()

    def __iter__(self):
        return self.read_chunks()

    def read(self):
        self._remove_initial_comments()
        header_data = self._buffer_type.read_header(self._file_obj)
        chunk = self._file_obj.read()
        chunk, _  = self.__add_newline_to_end(chunk, chunk.size)
        return self._buffer_type.from_raw_buffer(chunk)

    def read_chunks(self):
        self._remove_initial_comments()
        header_data = self._buffer_type.read_header(self._file_obj)
        chunk = self.__get_buffer()
        total_bytes = 0

        while not self._is_finished:
            total_bytes += chunk.size
            logger.debug(
                f"Read chunk of size {repr_bytes(chunk.size)} from {self._f_name}. (Total {repr_bytes(total_bytes)})"
            )
            buff = self._buffer_type.from_raw_buffer(chunk, header_data=header_data)
            self._file_obj.seek(buff.size - self._chunk_size, 1)
            yield wrapper(buff)
            chunk = self.__get_buffer()
        if chunk is not None and chunk.size:
            yield self._buffer_type.from_raw_buffer(chunk, header_data=header_data)

    def __add_newline_to_end(self, chunk, bytes_read):
        if chunk[bytes_read - 1] != ord("\n"):
            chunk = np.append(chunk, ord("\n"))
            bytes_read += 1
        if hasattr(self._buffer_type, "_new_entry_marker"):
            chunk = np.append(chunk, self._buffer_type._new_entry_marker)
            bytes_read += 1
        return chunk, bytes_read

    def __get_buffer(self):
        a, bytes_read = self.__read_raw_chunk()
        self._is_finished = bytes_read < self._chunk_size
        if bytes_read == 0:
            return None

        # Ensure that the last entry ends with newline. Makes logic easier later
        a, bytes_read = self.__add_newline_to_end(a, bytes_read)
        # if self._is_finished:
        #     if a[bytes_read - 1] != ord("\n"):
        #         a = np.append(a, ord("\n"))
        #         bytes_read += 1
        #     if hasattr(self._buffer_type, "_new_entry_marker"):
        #         a = np.append(a, self._buffer_type._new_entry_marker)
        #         bytes_read += 1
        return a[:bytes_read]

    def __read_raw_chunk(self):
        b = np.frombuffer(self._file_obj.read(self._chunk_size), dtype="uint8")
        return b, b.size
        # array = my_empty(self._chunk_size, dtype="uint8")
        # bytes_read = self._file_obj.readinto(array)
        # return array, bytes_read

    def _remove_initial_comments(self):
        if self._buffer_type.COMMENT == 0:
            return
        for line in self._file_obj:
            print(line)
            if line[0] != self._buffer_type.COMMENT:
                self._file_obj.seek(-len(line), 1)
                break


class NpBufferedWriter:
    """ 
    File writer that can write @npdataclass objects
    to file
    """

    def __init__(self, file_obj, buffer_type):
        self._file_obj = file_obj
        self._buffer_type = buffer_type
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file_obj:
            self._file_obj.close()

    def close(self):
        self._file_obj.close()

    def write(self, data: npdataclass):
        """Write the provided data to file

        Parameters
        ----------
        data : npdataclass
            Data set containing entries

        """
        if isinstance(data, NpDataclassStream):
            for buf in data:
                self.write(buf)
            return 
        bytes_array = self._buffer_type.from_data(data)
        self._file_obj.write(bytes(bytes_array))  # .tofile(self._file_obj)
        self._file_obj.flush()
        logger.debug(
            f"Wrote chunk of size {repr_bytes(bytes_array.size)} to {self._f_name}"
        )


def chunk_lines(stream, n_lines):
    cur_buffers = []
    remaining_lines = n_lines
    for chunk in stream:
        n_lines_in_chunk = len(chunk)
        while n_lines_in_chunk >= remaining_lines:
            cur_buffers.append(chunk[:remaining_lines])
            yield np.concatenate(cur_buffers)
            cur_buffers = []
            chunk = chunk[remaining_lines:]
            remaining_lines = n_lines
            n_lines_in_chunk = len(chunk)
        cur_buffers.append(chunk)
        remaining_lines -= n_lines_in_chunk
    yield np.concatenate(cur_buffers)
