import dataclasses
import pytest
import numpy as np
from bionumpy import AminoAcidEncoding, DNAEncoding, EncodedArray, BaseEncoding
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.bnpdataclass.bnpdataclass import make_dataclass, BNPDataClass, dynamic_concatenate
from bionumpy.bnpdataclass.bnpdataclassfunction import bnpdataclassfunction
from numpy.testing import assert_equal

from bionumpy.datatypes import SequenceID
from bionumpy.encodings.bool_encoding import bool_string
from bionumpy.util.testing import assert_bnpdataclass_equal
# import pandas as pd
import bionumpy as bnp


@bnpdataclass
class Person:
    name: str
    age: int


@pytest.fixture
def data():
    return Person(['knut', 'per', 'jon', 'erling'],
                  [10, 20, 30, 40])


@pytest.fixture
def data_list():
    return [Person(['knut', 'per', 'jon', 'erling'],
                   [10, 20, 30, 40])
            for _ in range(100)]


def test_add_fields():
    @bnpdataclass
    class BaseDC:
        sequence_aa: AminoAcidEncoding

    base_obj = BaseDC(['ACD', "EEA"])
    for field_map_dict in [{"sequence": DNAEncoding}, None]:
        res_obj = base_obj.add_fields({"sequence": ['AA', 'ACT']}, field_map_dict)

        print(res_obj)

        assert all(field.name in ['sequence', 'sequence_aa'] for field in dataclasses.fields(res_obj))
        assert res_obj.sequence.tolist() == ["AA", "ACT"]  # TODO: fix type hinting for fully dynamic stuff


def test_extend():
    @bnpdataclass
    class BaseDC:
        sequence_aa: AminoAcidEncoding

    extended_class = BaseDC.extend((('sequence', DNAEncoding), ('s1', int)))
    assert issubclass(extended_class, BaseDC)
    assert extended_class.__name__ == "DynamicBaseDC"
    assert all(field.name in ['sequence', 'sequence_aa', 's1'] for field in dataclasses.fields(extended_class))


def test_make_dataclass():
    new_cls = make_dataclass([("sequence", DNAEncoding), ('signal1', int)])

    assert issubclass(new_cls, BNPDataClass)
    assert new_cls.__name__ == "DynamicDC"
    assert all(field.name in ['sequence', 'signal1'] for field in dataclasses.fields(new_cls))


def add(a, b):
    return a + b


@pytest.mark.skip
def test_from_pandas():
    persons = {'name': ['knut', 'marit'],
               'age': [35, 30]}
    df = pd.DataFrame(persons)
    obj = Person(df.name, df.age)
    assert_bnpdataclass_equal(
        obj,
        Person(persons['name'], persons['age']))


@bnpdataclass
class MyClass:
    a: int
    b: int


def test_keyword_init():
    MyClass(a=[10, 20], b=[100, 200])


@pytest.mark.skip("not implemented")
def test_bnpdataclassfunction():
    bnp_add = bnpdataclassfunction("a", "b", (add))
    assert_equal(bnp_add(MyClass([10], [20])), [30])

# @pytest.mark.skip("Deprecated")
def test_set_get_context():
    data = MyClass(a=[10, 20], b=[100, 200])
    context = "Test test"
    data.set_context("test", context)
    assert data.get_context("test") == context


@pytest.mark.parametrize("file", [
    "variants.vcf",
    "variants_with_header.vcf"
])
def test_read_header(file,data_path):
    file = data_path/file
    chunks = list(bnp.open(file).read_chunks())
    true_header = "".join(line for line in open(file) if line.startswith("#"))
    for chunk in chunks:
        header = chunk.get_context("header")
        assert header == true_header


def test_dynamic_join(data_list):
    truth = np.concatenate(data_list)
    solution = dynamic_concatenate(data_list)
    assert_bnpdataclass_equal(truth, solution)


@bnpdataclass
class TwoPerson:
    person_1: Person
    person_2: Person
    relation: str


def test_hierachical(data):
    two_persons = TwoPerson(data, data, ['b', 'b', 'f', 'f'])
    assert_bnpdataclass_equal(two_persons.person_1, data)
    print(two_persons)


def test_tolist(data):
    entries = data.tolist()
    for entry in entries:
        assert isinstance(entry, Person.dataclass)
        assert isinstance(entry.name, str)
        assert isinstance(entry.age, int), type(entry.age)

@pytest.fixture()
def bool_class():
    @bnpdataclass
    class BNPDC:
        sequence_id: SequenceID
        test_field: bool_string

    return BNPDC

def test_bool_class(bool_class):
    obj = bool_class(sequence_id=['hei', 'ja'],
                     test_field=['True', 'False'])
    from bionumpy.io.delimited_buffers import DelimitedBuffer
    buffer = DelimitedBuffer.from_data(obj)



