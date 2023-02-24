import logging
import numpy as np
from typing import List, Union, Iterable, Tuple, Dict, Any
from ..bnpdataclass import replace
from .global_offset import GlobalOffset
from ..streams import groupby
from ..encodings.string_encodings import StringEncoding
from .genome_context_base import GenomeContextBase

logger = logging.getLogger(__name__)


class GenomeError(Exception):
    pass


class GenomeContext(GenomeContextBase):
    def __init__(self, chrom_size_dict: Dict[str, int], ignored=None):
        self._ignored = ignored
        if ignored is None:
            self._ignored = set([])
        keys = [name for name in chrom_size_dict if name not in ignored] + list(ignored)
        chrom_size_dict = {key: chrom_size_dict[key] for key in keys}
        self._included = [chrom for chrom in chrom_size_dict if chrom not in ignored]
        self._included_mask = np.array([chrom in self._included for chrom in chrom_size_dict])
        self._string_endcoding = StringEncoding(list(chrom_size_dict.keys()))
        self._chrom_size_dict = {key: value for key, value in chrom_size_dict.items() if key in self._included}
        self._global_offset = GlobalOffset(self._chrom_size_dict, string_encoding=self._string_endcoding)

    def __repr__(self):
        return repr(list(self._included)[:10] + ['...']*(len(self._included) > 10))

    @property
    def size(self):
        return sum(self._chrom_size_dict.values())

    @property
    def encoding(self):
        return self._string_endcoding

    @property
    def chrom_sizes(self):
        return self._chrom_size_dict

    @property
    def global_offset(self):
        return self._global_offset

    def is_included(self, chromosomes):
        return self._included_mask[chromosomes.raw()]

    def mask_data(self, data, chromosome_field_name='chromosome'):
        encoded_chromosomes = self.encoding.encode(getattr(data, chromosome_field_name))
        data = replace(data, **{chromosome_field_name: encoded_chromosomes})
        mask = self.is_included(encoded_chromosomes)
        return data[mask]

    @classmethod
    def from_dict(cls, chrom_size_dict):
        f = lambda key: '_' not in key
        return cls(chrom_size_dict, # {key: value for key, value in chrom_size_dict.items() if f(key)},
                   {key for key in chrom_size_dict if not f(key)})

    def chromosome_order(self):
        return (key for key in self._chrom_size_dict if '_' not in key)

    def is_compatible(self, other):
        return (self._chrom_size_dict == other._chrom_size_dict) and (self._included == other._included)

    def _included_groups(self, grouped):
        for name, group in grouped:
            if name in self._ignored:
                continue
            if name not in self._included:
                raise GenomeError(f'{name} not included in genome: {set(self._chrom_size_dict.keys())}')
            yield name, group

    def iter_chromosomes(self, data, dataclass, group_field='chromosome'):
        real_order = self.chromosome_order()
        grouped = groupby(data, group_field)
        grouped = self._included_groups(grouped)
        next_name, next_group = next(grouped, (None, None))
        seen = []
        seen_group = []
        for name in real_order:
            if name == next_name:
                seen_group.append(next_name)
                logger.debug(f'Yielding data for {name}')
                yield next_group
                next_name, next_group = next(grouped, (None, None))
                if next_name in seen:
                    raise GenomeError(
                        f'Sort order discrepancy ({next_name}): Genome so far {seen}, data: {seen_group}')
            else:
                logger.debug(f'Yielding empty data for {name}')
                yield dataclass.empty()
            seen.append(name)
        if next(grouped, None) is not None:
            raise GenomeError()
