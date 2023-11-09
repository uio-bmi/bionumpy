import numpy as np


class StringArray(np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, data):
        self._data = np.asanyarray(data, dtype='S')

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

    def _convert_input(self, value):
        if isinstance(value, str):
            return bytes(value)
        elif isinstance(value, self.__class__):
            return value.raw()
        return np.asanyarray(value, dtype='S')
