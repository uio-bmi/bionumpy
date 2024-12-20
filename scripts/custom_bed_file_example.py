import bionumpy as bnp
from bionumpy import BedBuffer

'''
This script shows how to create a custom BED format for reading a BED file with an additional sequence column.
The BED file is assumed to have the following format:
    - chromosome
    - start
    - end
    - sequence

By creating a custom buffer class, we can specify the dataclass that should be used to store the data.
'''
filepath = "example_data/interval_with_sequence.bed"


@bnp.bnpdataclass.bnpdataclass
class BedWithSequence(bnp.datatypes.Interval):
    sequence: str


class CustomBedBuffer(BedBuffer):
    dataclass = BedWithSequence


def test_bed():
    # Reading the file with the normal bed buffer looses the sequence field:
    bed_wo_sequence = bnp.open(filepath).read()
    print(bed_wo_sequence)
    '''
    Interval with 2 entries
               chromosome                    start                     stop
                     chr1                        0                       10
                     chr1                        0                       20
    '''
    bed = bnp.open(filepath, buffer_type=CustomBedBuffer).read()
    print(bed)
    '''
    BedWithSequence with 2 entries
               chromosome                    start                     stop                 sequence
                     chr1                        0                       10                     GCCT
                     chr1                        0                       20                     ACGT
'''
