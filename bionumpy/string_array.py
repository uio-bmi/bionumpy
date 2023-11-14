import numpy as np

class StringArray(np.lib.mixins.NDArrayOperatorsMixin):
    wrapped_functions = ['size', 'shape', 'ndim', '__len__']

    def __init__(self, data):
        self._data = np.asanyarray(data, dtype='S')

    def __repr__(self):
        if self._data.ndim>=1:
            return repr(self._data[:10].astype('U'))
        else:
            return self._data.astype('U').__repr__()

    def raw(self):
        return self._data

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handle numpy ufuncs called on EnocedArray objects

        Only support euqality checks for now. Numeric operations must be performed
        directly on the underlying data.
        """

        if method == "__call__" and ufunc.__name__ in ("equal", "not_equal"):
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

    def __setitem__(self, item, value):
        self._data[item] == self._convert_input(item, value)

    def __getattr__(self, name):
        if name in self.wrapped_functions:
            return getattr(self._data, name)
        return getattr(self, name)

    def _convert_input(self, value):
        if isinstance(value, str):
            return bytes(value)
        elif isinstance(value, self.__class__):
            return value.raw()
        return np.asanyarray(value, dtype='S')

    def __len__(self):
        return len(self._data)


def string_array(input_data):
    print(type(input_data))
    if isinstance(input_data, (list, np.ndarray)):
        return StringArray(input_data)
    elif isinstance(input_data, StringArray):
        return input_data.copy()
    else:
        raise TypeError(f'Cannot convert {input_data} to StringArray')


def as_string_array(input_data):
    if isinstance(input_data, StringArray):
        return input_data
    return string_array(input_data)