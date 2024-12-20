import dataclasses
import inspect
import logging
from collections import defaultdict
from typing import List, Type, Dict, Iterable, Union, Optional, Any, Tuple
from numpy.typing import ArrayLike
from npstructures.npdataclasses import npdataclass, NpDataClass, shallow_tuple
from npstructures import RaggedArray
import numpy as np

from ..typing import SequenceID
from .pandas_adaptor import pandas_adaptor
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encoded_array import as_encoded_array
from ..encodings import Encoding, NumericEncoding
from ..encodings.alphabet_encoding import FlatAlphabetEncoding
from ..string_array import as_string_array
from ..util import is_subclass_or_instance

logger = logging.getLogger(__name__)


def get_vanilla_generator(object):
    if isinstance(object, np.ndarray):
        convertor = (lambda x: x.item()) if (object.ndim == 1) else (lambda x: x.tolist())
        return (convertor(o) for o in object)
    if isinstance(object, (EncodedArray, EncodedRaggedArray)):
        return (c.to_string() for c in object)
    if isinstance(object, RaggedArray):
        return (row.tolist() for row in object)
    if isinstance(object, BNPDataClass):
        return object.toiter()

    return (c.to_string() for c in object)


class BNPDataClass(NpDataClass):

    def todict(self) -> Dict[str, ArrayLike]:
        '''
        Convert the data into a dictionary with the field names as keys and the corresponding data as values.

        Returns
        -------
            dict[str, ArrayLike]
        '''
        field_dict = {}
        for field in dataclasses.fields(self):
            pandas_obj = pandas_adaptor.pandas_converter(getattr(self, field.name))
            if isinstance(pandas_obj, dict):
                field_dict.update({f'{field.name}.{k}': v for k, v in pandas_obj.items()})
            else:
                field_dict[field.name] = pandas_obj

        return field_dict

    def topandas(self) -> 'pandas.DataFrame':
        '''
        Convert the data into a pandas DataFrame with the fields as columns

        Returns
        -------
            pandas.DataFrame
        '''
        return pandas_adaptor.get_data_frame(self.todict())
        # return pd.DataFrame(self.todict())

    @classmethod
    def from_data_frame(cls, df: 'pandas.DataFrame') -> 'BNPDataClass':
        '''
        Convert a pandas DataFrame into a BNPDataClass object.
        The columns of the dataframe are used as fields in the BNPDataClass object.

        Parameters
        ----------
        df: pandas.DataFrame

        Returns
        -------
            BNPDataClass
        '''
        d = df.to_dict('series')
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, dict_object: Dict[str, Any]) -> 'BNPDataClass':
        '''
        Convert a dictionary into a BNPDataClass object. The keys of the dictionary are used as field names in the BNPDataClass object.
        Parameters
        ----------
        dict_object: dict

        Returns
        -------
            BNPDataClass
        '''
        dict_names = [name.split('.')[0] for name in dict_object.keys()]
        field_names = {field.name for field in dataclasses.fields(cls)}
        logger.info(f'Dropping columns: {[n for n in dict_names if n not in field_names]}')
        new_dict = defaultdict(dict)
        for name, value in dict_object.items():
            if '.' in name:
                name, subname = name.split('.', maxsplit=1)
                new_dict[name][subname] = value
            else:
                new_dict[name] = value
        for field in dataclasses.fields(cls):
            if isinstance(new_dict[field.name], dict):
                assert is_subclass_or_instance(field.type, BNPDataClass), field
                new_dict[field.name] = field.type.from_dict(new_dict[field.name])
        return cls(**new_dict)

    def tolist(self)-> List['dataclass']:
        """
        Convert the data into a list of entries from the
        corrsponding dataclass with normal python types.
        Similar to np.tolist and pd.tolist.
        This is good for debugging, but for real applications
        requires a lot of memory allocation. For iterating over
        the data, use `toiter` instead.

        Returns
        -------
            list[cls.dataclass]
        """
        return list(self.toiter())
        lists = tuple(f.tolist() for f in shallow_tuple(self))
        return list(self.dataclass(*row) for row in zip(*lists))

    def toiter(self) -> Iterable['dataclass']:
        """
        Convert the data into an iterator of entries from the
        corrsponding dataclass with normal python types.

        Returns
        -------
            Iterable[cls.dataclass]

        """
        iters = tuple(get_vanilla_generator(f)
                      for f in shallow_tuple(self))
        return (self.dataclass(*row) for row in zip(*iters))

    to_iter = toiter

    @classmethod
    def extend(cls, fields: tuple, name: str = None) -> Type['BNPDataClass']:
        """
        Parameters
        ----------
        fields: tuple
            A tuple in format (field_name, field_type) for the new fields to be added
        name: str
            The optional user-defined name for the new class

        Returns
        --------
        BNPDataClass with added fields

        Examples
        ---------

        >>> from bionumpy.bnpdataclass import bnpdataclass
        >>> from bionumpy.encodings import AminoAcidEncoding, DNAEncoding
        >>> @bnpdataclass
        ... class BaseDC:
        ...     sequence_aa: AminoAcidEncoding

        >>> extended_class = BaseDC.extend((('sequence', DNAEncoding), ('s1', int)))
        >>> assert all(field.name in ['sequence', 'sequence_aa', 's1'] for field in dataclasses.fields(extended_class))
        >>> print([field.name for field in dataclasses.fields(extended_class)])
        ['sequence_aa', 'sequence', 's1']

        """
        cls_name = name if name is not None else f"Dynamic{cls.__name__}" if cls.__name__[
                                                                             :7] != 'Dynamic' else cls.__name__
        return bnpdataclass(dataclasses.make_dataclass(cls_name, bases=(cls,), fields=fields))

    def add_fields(self, fields: Dict[str, ArrayLike], field_type_map: dict = None) -> 'BNPDataClass':
        """
        Parameters
        ----------
        fields: dict
            a dictionary in containing the names of the new fields as keys and lists of values for each of the field as values

        field_type_map: dict
            a dictionary with field names as keys and types as values; for basic types, they can be inferred from the data and don't need to be
            specified; but for fields that need to use some of the encodings, the specific encoding can be provided here

        Returns
        --------
        BNPDataClass object with added fields with the provided values

        Examples
        ---------

        >>> from bionumpy.bnpdataclass import bnpdataclass
        >>> from bionumpy.encodings import AminoAcidEncoding, DNAEncoding
        >>> @bnpdataclass
        ... class BaseDC:
        ...     sequence_aa: AminoAcidEncoding

        >>> base_obj = BaseDC(['ACD', "EEA"])
        >>> res_obj = base_obj.add_fields({"sequence": ['AA', 'ACT']}, field_type_map={'sequence': DNAEncoding})
        >>> print(res_obj)
        DynamicBaseDC with 2 entries
                      sequence_aa                 sequence
                              ACD                       AA
                              EEA                      ACT

        """
        for name in fields.keys():
            if not name.isidentifier():
                raise TypeError(f"Field name must be a valid identifier (No whitespace or dots and such): {name}")
        fields_with_types = _extract_field_types(fields, field_type_map)
        new_class = self.__class__.extend(tuple(fields_with_types.items()))
        return new_class(**{**vars(self), **fields})

    @classmethod
    def from_entry_tuples(cls, tuples: Iterable[tuple]) -> 'BNPDataClass':
        return cls(*(list(c) for c in zip(*tuples)))

    def sort_by(self, field_name: str) -> 'BNPDataClass':
        """
        Sort the data by the given field

        Parameters
        ----------
        field_name: str

        Returns
        -------
        BNPDataClass

        """
        return self[np.argsort(getattr(self, field_name))]

    def set_context(self, name: str, value: Any):
        """
        Set a context value for the object, typycally used for storing auxillary information like header information
        Parameters
        ----------
        name: str
        value: Any

        Returns
        -------

        """
        if not hasattr(self, '_context'):
            self._context = dict()
        self._context[name] = value

    def get_context(self, name: str)->Any:
        """
        Get a context value for the object, typycally used for storing auxillary information like header information
        Parameters
        ----------
        name: str

        Returns
        -------

        """
        logger.warning(f'Deprecated method set_context in BNPDataClass')
        if not hasattr(self, '_context'):
            self._context = dict()
        return self._context[name]

    def has_context(self, name):
        return hasattr(self, '_context') and name in self._context


def bnpdataclass(base_class: type) -> Type[BNPDataClass]:
    """Create a `bnpdataclass` from a class with fields specified

    A wrapper around `@npdataclass` that includes implicit format conversion
    for strings and other encoded data. `@npdataclass` is again a wrapper
    around `dataclasses.dataclass` but where all the fields are assumed to be
    objects that supports advanced numpy indexing so that any dataclass objects
    are also indexible like a numpy array.

    `bnpdataclass` classes are meant to be dataclasses and so should not have 
    any methods other than those implicitly given by `npdataclass`.

    Parameters
    ----------
    base_class : type
        Base class that defines the fields of the dataclass.

    Returns
    -------
    bnpdataclass
        `bnpdataclass` object that supports numpy like indexing

    Examples
    --------
    >>> from bionumpy.bnpdataclass import bnpdataclass
    >>> @bnpdataclass
    ... class Person:
    ...       name: str
    ...       age: int
    ... 
    >>> data = Person(["Knut", "Ivar", "Geir"], [35, 30, 40])
    >>> print(data)
    Person with 3 entries
                         name                      age
                         Knut                       35
                         Ivar                       30
                         Geir                       40

    >>> print(data[[0,2]])
    Person with 2 entries
                         name                      age
                         Knut                       35
                         Geir                       40
    
    >>> print(data[[False,True, False]])
    Person with 1 entries
                         name                      age
                         Ivar                       30
    """

    class NewClass(npdataclass(base_class), BNPDataClass):
        @classmethod
        def _implicit_format_conversion(cls, obj: npdataclass):
            """Convert the data in given in the init into numpy like data

            This is convenience functionionality that converts e.g. list of strings
            into EncodedRaggedArray objects, and numeric data into `np.ndarray` objects.

            Called by the `__init__` method from `npdataclass`.

            Parameters
            ----------
            cls : 3
                The class this is called from
            obj : npdataclass
                The partially initialzed object from `npdataclass` `__init__`

            """
            for field in dataclasses.fields(obj):
                pre_val = getattr(obj, field.name)
                try:
                    val = cls.__convert_single_field(field, pre_val)
                except Exception as e:
                    raise ValueError(f"Error when converting {field.name} to {field.type} with value {pre_val}") from e

                setattr(obj, field.name, val)

        @classmethod
        def __convert_single_field(cls, field, pre_val):
            numeric_types = (int, float, bool)
            optional_numeric_types = tuple(Optional[t] for t in numeric_types)
            if field.type == Union[BNPDataClass, str]:
                if isinstance(pre_val,
                              (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)) or \
                        hasattr(pre_val, 'to_numpy'):
                    val = as_encoded_array(pre_val)
                elif True or isinstance(pre_val, BNPDataClass):
                    val = pre_val
                else:
                    assert False, (field.type, type(pre_val))

            elif field.type in numeric_types + optional_numeric_types:
                val = np.asanyarray(pre_val)
            elif field.type == str:
                assert isinstance(pre_val, (
                    str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)) or hasattr(pre_val,
                                                                                                      'to_numpy'), (
                    field, pre_val, type(pre_val))
                val = as_encoded_array(pre_val)
            elif field.type == SequenceID or field.type == List[str]:
                if isinstance(pre_val, EncodedArray):
                    val = pre_val
                else:
                    val = as_string_array(pre_val)
            elif is_subclass_or_instance(field.type, Encoding):
                if is_subclass_or_instance(field.type, NumericEncoding):
                    assert isinstance(pre_val,
                                      (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)), \
                        (field, pre_val, type(pre_val))
                    val = as_encoded_array(pre_val, field.type)
                elif getattr(field.type, 'returns_raw', False) and isinstance(pre_val, (np.ndarray, np.generic)):
                    val = pre_val
                else:
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, bool)) or hasattr(pre_val,
                                                                                                               'to_numpy'), (
                    field, pre_val, type(pre_val), isinstance(pre_val, np.generic))
                    val = as_encoded_array(pre_val, field.type)
                # must do as_encoded and not explicit encode as pre_val might already
                # be encoded
                if isinstance(field.type, FlatAlphabetEncoding):
                    val = val.ravel()
            elif field.type == List[int] or field.type == List[bool] or field.type == List[float]:
                if not isinstance(pre_val, RaggedArray):
                    try:
                        val = RaggedArray(pre_val)
                    except TypeError as e:
                        val = np.asanyarray(pre_val)
                else:
                    val = pre_val
            elif inspect.isclass(field.type) and issubclass(field.type, BNPDataClass):
                # assert isinstance(pre_val, (field.type, field.type._single_entry)), (field.type, type(pre_val))
                val = pre_val
            else:
                assert False, field.type
            return val

    NewClass.__name__ = base_class.__name__
    NewClass.__qualname__ = base_class.__qualname__
    NewClass.__doc__ = dataclasses.dataclass(base_class).__doc__
    return NewClass


def make_dataclass(fields: List[Tuple], name: str = "DynamicDC", bases=()) -> Type[BNPDataClass]:
    """
    Constructs a dynamic dataclass from a list of attributes

    Parameters
    ----------
    fields: list
        a list of tuples in format (field_name, field_type) to be used to construct the dynamic bnp dataclass

    name: str
        optional name of new class

    Returns
    -------
    new BNPDataClass

    """
    return bnpdataclass(dataclasses.make_dataclass(name, fields=fields, bases=bases))


def narrow_type(bnp_dc: Type[BNPDataClass], field_name: str, field_type: type):
    """
    Resticts the type of a field in a BNPDataClass

    Parameters
    ----------
    bnp_dc: Type[BNPDataClass]
    field_name: str
    field_type: type

    Returns
    -------
    Type[BNPDataClass]

    """
    new_fields = [(f.name, field_type) if f.name==field_name else (f.name, f.type, f) for f in dataclasses.fields(bnp_dc)]
    return make_dataclass(new_fields, name=bnp_dc.__name__, bases=(bnp_dc,))


def _extract_field_types(fields_with_values: dict, field_type_map: dict = None) -> dict:
    fields = {}
    for field_name in fields_with_values.keys():
        _assert_all_same_type(fields_with_values[field_name])

        if field_type_map is not None and field_name in field_type_map:
            field_type = field_type_map[field_name]
        elif isinstance(fields_with_values[field_name][0], EncodedArray):
            field_type = type(fields_with_values[field_name][0].encoding)
        else:
            field_type = type(fields_with_values[field_name][0])
        if fields_with_values[field_name] is not None:
            fields[field_name] = field_type

    return fields


def _assert_all_same_type(values):
    original_type = type(values[0])
    assert all(isinstance(val, original_type) for val in values), (original_type, [type(val) for val in values])


def dynamic_concatenate(dataclass_iter: Iterable[BNPDataClass]):
    iterable = iter(dataclass_iter)

    first = next(iterable)
    first_class = first.__class__
    fields = [[vals] for vals in first.shallow_tuple()]
    l = len(first)
    for c in iterable:
        for f, vals in zip(fields, c.shallow_tuple()):
            f.append(vals)
        l += len(c)
        print(l)
    print('Joining fields', sum(len(f) for f in fields[0]))
    for i, f in enumerate(fields):
        fields[i] = np.concatenate(f)
    print('creating object')
    return first_class(*fields)

