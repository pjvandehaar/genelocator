
import bisect
import intervaltree

class GeneLocator:
    def __init__(self, genes):
        '''genes is like [{chrom:"1", start:123, end:234, ensg:"ENSG00345", symbol:"ACG4"},...]'''
        self._gene_info = {}
        self._its = {} # an interval tree for each chromosome
        self._gene_starts = {} # a list of (start, ensg) pairs for each chromosome
        self._gene_ends = {} # a list of (end, ensg) pairs for each chromosome
        for gene in genes:
            chrom = gene['chrom'][len('chr'):] if gene['chrom'].startswith('chr') else gene['chrom']
            start = gene['start']; end = gene['end']; ensg = gene['ensg']; symbol = gene['symbol']
            if not isinstance(start, int) or not isinstance(end, int): raise Exception("start and end must be int, unlike {!r} and {!r}".format(start, end))
            if ensg in self._gene_info: raise Exception("The gene {!r} appears multiple times in this genes list (as {!r} and {!r})".format(ensg, self._gene_info[ensg], gene))
            self._gene_info[ensg] = (chrom, start, end, symbol)
            self._its.setdefault(chrom, intervaltree.IntervalTree()).addi(start, end, ensg)
            self._gene_starts.setdefault(chrom, []).append((start, ensg))
            self._gene_ends.setdefault(chrom, []).append((end, ensg))
        for chrom in self._its:
            self._gene_starts[chrom] = BisectFinder(self._gene_starts[chrom])
            self._gene_ends[chrom] = BisectFinder(self._gene_ends[chrom])

    def at(self, chrom, pos):
        if chrom.startswith('chr'): chrom = chrom[len('chr'):]
        if chrom == 'MT': chrom = 'M'
        if chrom not in self._its: raise Exception("Unknown chromosome: {!r}".format(chrom))
        # If any genes overlap this position, return them all.
        overlapping_genes = self._its[chrom].at(pos) if hasattr(self._its[chrom], 'at') else self._its[chrom].search(pos) # support intervaltree 2.x and 3.x
        if overlapping_genes: return [self._with_info(g.data) for g in overlapping_genes]
        # If only one direction has genes (either before or after this position), return the gene in that direction.
        prev_gene_end = self._gene_ends[chrom].get_item_before_or_at(pos)
        next_gene_start = self._gene_starts[chrom].get_item_after_or_at(pos)
        if prev_gene_end is None and next_gene_start is None: raise Exception('The position chr{!r}:{!r} has no genes before it or after it'.format(chrom, pos))
        if next_gene_start is None: return [self._with_info(prev_gene_end[1])]
        if prev_gene_end is None: return [self._with_info(next_gene_start[1])]
        # If both directions have genes, return whichever is closer.
        dist_to_prev_gene_end = abs(prev_gene_end[0] - pos)
        dist_to_next_gene_start = abs(next_gene_start[0] - pos)
        return [self._with_info(prev_gene_end[1])] if dist_to_prev_gene_end < dist_to_next_gene_start else [self._with_info(next_gene_start[1])]

    def _with_info(self, ensg):
        info = self._gene_info[ensg]
        return {'chrom':info[0], 'start':info[1], 'end':info[2], 'ensg':ensg, 'symbol': info[3]}

class BisectFinder(object):
    '''Given a list like [(123, 'foo'), (125, 'bar')...], BisectFinder helps you find the things before (or at) and after (or at) 124.'''
    def __init__(self, tuples):
        '''tuples is like [(123, 'foo'),...]'''
        tuples = sorted(tuples, key=lambda t:t[0])
        self._nums, self._values = list(zip(*tuples))
    def get_item_before_or_at(self, pos):
        '''If we get an exact match, let's return it'''
        idx = bisect.bisect_right(self._nums, pos) - 1 # note: bisect_{left,right} just deals with ties.
        if idx < 0: return None # It's fallen off the beginning!
        return (self._nums[idx], self._values[idx])
    def get_item_after_or_at(self, pos):
        if pos > self._nums[-1]: return None # it's fallen off the end!
        idx = bisect.bisect_left(self._nums, pos)
        return (self._nums[idx], self._values[idx])
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_before_or_at(22) is None
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_before_or_at(23) == (23,'foo')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_before_or_at(24) == (23,'foo')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_before_or_at(25) == (25,'bar')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_before_or_at(26) == (25,'bar')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_after_or_at(22) == (23,'foo')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_after_or_at(23) == (23,'foo')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_after_or_at(24) == (25,'bar')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_after_or_at(25) == (25,'bar')
assert BisectFinder([(23,'foo'),(25,'bar')]).get_item_after_or_at(26) is None
