import numpy as np
from numpy.typing import ArrayLike
from typing import List, Dict, Tuple, Optional
from .strops import ints_to_strings, int_lists_to_strings, float_to_strings
from ..encoded_array import EncodedArray, Encoding, EncodedRaggedArray, BaseEncoding, encoded_array_from_nparray, \
    change_encoding, as_encoded_array
from ..encodings.bool_encoding import bool_string
from ..encodings.string_encodings import StringEncoding
from npstructures import RaggedArray

from ..typing import SequenceID
from ..util import is_subclass_or_instance


def str_func(column):
    if column.encoding == BaseEncoding:
        return column
    elif isinstance(column.encoding, StringEncoding):
        return column.encoding.decode(column)
    else:
        return change_encoding(column, BaseEncoding)
    assert False, column.encoding


def str_matrix_func(column):
    n_rows, n_cols = column.shape
    a = column.as_bytes().reshape(n_rows * n_cols, -1)
    tabs = np.full((n_rows * n_cols, 1), ord("\t"))
    b = np.hstack([a, tabs])
    b = b.reshape((n_rows, -1))[:, :-1]
    return EncodedRaggedArray(EncodedArray(b.ravel(), BaseEncoding), np.full(b.shape[0], b.shape[-1]))


def seq_id_func(column):
    if isinstance(column, EncodedArray):
        if isinstance(column.encoding, StringEncoding):
            return column.encoding.decode(column)
    return encoded_array_from_nparray(column)


def optional_ints_to_strings(number: np.ndarray, missing_string='.') -> EncodedRaggedArray:
    if np.all(number) == np.nan:
        return as_encoded_array([missing_string] * len(number))
    return ints_to_strings(number)


def get_column(values, field_type) -> EncodedRaggedArray:
    def get_func_for_datatype(datatype):
        funcs = {int: ints_to_strings,
                 Optional[int]: optional_ints_to_strings,
                 str: str_func,  # lambda x: x,
                 bool_string: lambda x: bool_string.decode(x),
                 SequenceID: seq_id_func,
                 List[int]: int_lists_to_strings,
                 float: float_to_strings,
                 List[bool]: lambda x: int_lists_to_strings(x.astype(int), sep=""),
                 bool: lambda x: ints_to_strings(x.astype(int)),
                 List[str]: str_matrix_func
                 }
        if is_subclass_or_instance(datatype, Encoding) and not datatype==bool_string:
            encoding = datatype

            def dynamic(x):
                if isinstance(x, EncodedRaggedArray):
                    return EncodedRaggedArray(EncodedArray(encoding.decode(x.ravel()), BaseEncoding), x.shape)
                return EncodedArray(encoding.decode(x), BaseEncoding)

            return dynamic
        else:
            return funcs[datatype]

    return get_func_for_datatype(field_type)(values)


def dump_csv(data_dict: List[Tuple], sep: str = "\t") -> EncodedArray:
    """Put each field of the dataclass into a column in a buffer.

    Parameters
    data_dict: List[Tuple]
        A list of tuples where each tuple contains the field name and the field value.
    sep: str
        The separator to use between fields.

    Returns
    -------
    EncodedArray
        A buffer containing the data in CSV format.

    """

    columns = [get_column(value, key) for key, value in data_dict]
    lines = join_columns(columns, sep)
    return lines.ravel()


def join_columns(columns: List[EncodedRaggedArray], sep: str) -> EncodedRaggedArray:
    """
    Join columns into a single buffer.

    Parameters
    ----------
    columns: List[EncodedRaggedArray]
    sep: str

    Returns
    -------
    EncodedRaggedArray
        The lines of the buffer

    """
    lengths = np.concatenate([((column.lengths if
                                isinstance(column, RaggedArray)
                                else np.array([
                                                  column.shape[-1] if len(column.shape) == 2 else 1
                                              ] * len(column))) + 1
                               )[:, np.newaxis]

                              for column in columns], axis=-1).ravel()
    lines = EncodedRaggedArray(EncodedArray(np.empty(lengths.sum(), dtype=np.uint8), BaseEncoding),
                               lengths)
    n_columns = len(columns)
    for i, column in enumerate(columns):
        if (not isinstance(column, RaggedArray)) and len(column.shape) == 1:
            column = EncodedRaggedArray.from_numpy_array(column[:, np.newaxis])
        lines[i::n_columns, :-1] = column
    lines[:, -1] = sep
    lines[(n_columns - 1)::n_columns, -1] = "\n"
    return lines
