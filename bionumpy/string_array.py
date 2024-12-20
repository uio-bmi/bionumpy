import numpy as np

from .encoded_array import EncodedRaggedArray, BaseEncoding, EncodedArray, as_encoded_array, change_encoding


class StringArray(np.lib.mixins.NDArrayOperatorsMixin):
    """Wrapper around NumPy arrays of strings. Can be used as datatype in BNPDataClass fields."""
    wrapped_functions = ['size', 'shape', 'ndim', '__len__']
    wrapped_properies = ['T']

    def __init__(self, data):
        self._data = np.asanyarray(data, dtype='S')

    @property
    def encoding(self):
        return None

    def _as_bytes(self):
        data = self._data
        if not data.flags['C_CONTIGUOUS']:
            data = data.flatten()
        return data.view(np.uint8).reshape(data.shape + (-1,))

    def as_bytes(self):
        return self._as_bytes()

    def ravel(self):
        raveled = self._as_bytes().ravel()
        return self.__class__(raveled[raveled!=0].view('S1'))

    @property
    def lengths(self):
        return np.count_nonzero(self._as_bytes(), axis=-1)

    def __repr__(self):

        ndim = self._data.ndim
        if ndim == 1:
            return '\n'.join(b.decode() for b in self._data[:5].tolist())
        elif ndim == 2:
            return '\n'.join('\t'.join(b.decode() for b in line) for line in self._data[:5].tolist())
        else:
            return self._data.tolist().decode()

    def raw(self):
        return self._data

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handle numpy ufuncs called on EnocedArray objects
        Only support euqality checks for now. Numeric operations must be performed
        directly on the underlying data.
        """

        if method != "__call__":
            return NotImplemented
        if ufunc.__name__ in ("equal", "not_equal"):
            inputs = [self._convert_input(input) for input in inputs]
            assert len(inputs) == 2, 'only binary operations are supported'
            if ufunc.__name__ == 'equal':
                return inputs[0] == inputs[1]
            elif ufunc.__name__ == 'not_equal':
                return inputs[0] != inputs[1]
            return ufunc(*inputs)

        return NotImplemented

    def __getitem__(self, item):
        return self.__class__(self._data[item])

    def __array_function__(self, func, types, args, kwargs):
        if func == np.isin:
            if type(args[1]) == list and all(type(a) == str for a in args[1]):
                vals = as_string_array(args[1])
            elif isinstance(args[1], StringArray):
                vals = args[1]
            else:
                return NotImplemented
            return np.isin(self._data, vals._data)
        if not all(issubclass(t, StringArray) for t in types):
            return NotImplemented
        if func == np.concatenate:
            return self.__class__(func([a._data for a in args[0]]))
        return NotImplemented


    def __setitem__(self, item, value):
        self._data[item] == self._convert_input(item, value)

    def __getattr__(self, name):
        if name in self.wrapped_functions:
            return getattr(self._data, name)
        if name in self.wrapped_properies:
            return self.__class__(getattr(self._data, name))
        raise AttributeError(f'{self.__class__.__name__} has no attribute {name}')

    def _convert_input(self, value):
        if isinstance(value, str):
            return bytes(value, 'ascii')
        elif isinstance(value, self.__class__):
            return value.raw()
        elif isinstance(value, (EncodedArray, EncodedRaggedArray)):
            print(value)
            return string_array(value)
        return np.asanyarray(value, dtype='S')

    def tolist(self):
        byte_list = self._data.tolist()
        if isinstance(byte_list, bytes):
            return byte_list.decode()
        return [s.decode() for s in byte_list]

    to_string=tolist

    def __len__(self):
        return len(self._data)

def printio(func):
    def new_func(*args, **kwargs):
        print(args)
        res = func(*args, **kwargs)
        print('>', res)
        return res
    return new_func

# @printio
def string_array(input_data):
    if isinstance(input_data, list) and len(input_data) > 0 and isinstance(input_data[0], StringArray):
        return string_array([i.raw() for i in input_data])
    if isinstance(input_data, (list, str)):
        return StringArray(np.array(input_data, dtype='S'))
    elif isinstance(input_data, np.ndarray):
        return StringArray(input_data)
    elif isinstance(input_data, StringArray):
        return input_data.copy()
    elif isinstance(input_data, (EncodedRaggedArray, EncodedArray)):
        if not input_data.encoding == BaseEncoding:
            input_data = input_data.encoding.decode(input_data)
        array = input_data.raw()
        if isinstance(input_data, EncodedRaggedArray):
            if len(input_data) == 0:
                return StringArray(np.array([], dtype='S'))
            array = array.as_padded_matrix(side='right')
        n_bytes = array.shape[-1]
        return StringArray(array.flatten().view(f'|S{n_bytes}'))

    if hasattr(input_data, 'to_numpy'):
        return string_array(input_data.to_numpy().tolist())
    else:
        raise TypeError(f'Cannot convert {input_data} to StringArray ({type(input_data)})')


def as_string_array(input_data):
    if isinstance(input_data, StringArray):
        return input_data
    return string_array(input_data)