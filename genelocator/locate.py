"""Helpers for nearest gene location"""

import bisect
import re
import typing as ty

import intervaltree

from . import exception as gene_exc


CHROM_PREFIX = re.compile('^chr')


def _chrom_helper(value: str) -> str:
    """
    A lot of people add a prefix when referring to chromosomes; standardize chr1 and 1 to the latter.
    This makes it easier to run queries
    """
    return CHROM_PREFIX.sub('', value)


class GeneLocator:
    def __init__(self, genes: ty.Iterable[dict]):
        """genes is like [{chrom: "1", start: 123, end: 234, ensg: "ENSG00345", symbol: "ACG4"},...]"""
        self._gene_info = {}
        self._its = {}  # an interval tree for each chromosome
        self._gene_starts = {}  # a searchable list of (start, ensg) pairs, per chromosome
        self._gene_ends = {}  # a list of (end, ensg) pairs for each chromosome

        for gene in genes:
            chrom = _chrom_helper(gene['chrom'])
            start = gene['start']
            end = gene['end']
            ensg = gene['ensg']
            symbol = gene['symbol']
            if not isinstance(start, int) or not isinstance(end, int):
                raise gene_exc.LookupCreateError("start and end must be int, unlike {!r} and {!r}".format(start, end))
            if ensg in self._gene_info:
                raise gene_exc.LookupCreateError(
                    "The gene {!r} appears multiple times in this genes list (as {!r} and {!r})".format(
                        ensg, self._gene_info[ensg], gene))
            self._gene_info[ensg] = (chrom, start, end, symbol)
            self._its.setdefault(chrom, intervaltree.IntervalTree()).addi(start, end, ensg)
            self._gene_starts.setdefault(chrom, []).append((start, ensg))
            self._gene_ends.setdefault(chrom, []).append((end, ensg))

        for chrom in self._its:
            self._gene_starts[chrom] = BisectFinder(self._gene_starts[chrom])
            self._gene_ends[chrom] = BisectFinder(self._gene_ends[chrom])

    def at(self, chrom, pos, *, strict=True):
        """Locate a gene from position coordinates"""

        chrom = _chrom_helper(chrom)

        if chrom == 'MT':
            # FIXME: Should we be coercing chromosome names in a generic reusable class?
            chrom = 'M'

        if chrom not in self._its:
            if strict:
                raise gene_exc.BadCoordinateException("Unknown chromosome: {!r}".format(chrom))

        # If any genes overlap this position, return them all.
        overlapping_genes = self._its[chrom].at(pos)
        if overlapping_genes:
            return [self._serialize(g.data) for g in overlapping_genes]

        # If only one direction has genes (either before or after this position), return the gene in that direction.
        prev_gene_end = self._gene_ends[chrom].get_item_before_or_at(pos)
        next_gene_start = self._gene_starts[chrom].get_item_after_or_at(pos)
        if prev_gene_end is None and next_gene_start is None:
            raise gene_exc.NoResultsFoundException(
                'The position chr{!r}:{!r} has no genes before it or after it'.format(chrom, pos))
        if next_gene_start is None:
            return [self._serialize(prev_gene_end[1])]
        if prev_gene_end is None:
            return [self._serialize(next_gene_start[1])]

        # If both directions have genes, return whichever is closer.
        dist_to_prev_gene_end = abs(prev_gene_end[0] - pos)
        dist_to_next_gene_start = abs(next_gene_start[0] - pos)
        if dist_to_prev_gene_end < dist_to_next_gene_start:
            return [self._serialize(prev_gene_end[1])]
        else:
            return [self._serialize(next_gene_start[1])]

    def _serialize(self, ensg) -> dict:
        """Return a serialized representation of the gene data"""
        info = self._gene_info[ensg]
        return {'chrom': info[0], 'start': info[1], 'end': info[2], 'ensg': ensg, 'symbol': info[3]}


class BisectFinder:
    """
    Given a list like [(123, 'foo'), (125, 'bar')...], BisectFinder helps you find the things before (or at)
     and after (or at) 124.
    """
    def __init__(self, tuples):
        """tuples is like [(123, 'foo'),...]"""
        tuples = sorted(tuples, key=lambda t: t[0])
        self._nums, self._values = list(zip(*tuples))

    def get_item_before_or_at(self, pos):
        """If we get an exact match, let's return it"""
        idx = bisect.bisect_right(self._nums, pos) - 1  # note: bisect_{left,right} just deals with ties.
        if idx < 0:
            return None  # It's fallen off the beginning!
        return self._nums[idx], self._values[idx]

    def get_item_after_or_at(self, pos):
        if pos > self._nums[-1]:
            return None  # it's fallen off the end!
        idx = bisect.bisect_left(self._nums, pos)
        return self._nums[idx], self._values[idx]
