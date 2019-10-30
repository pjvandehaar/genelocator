import pytest

from genelocator import get_genelocator, exceptions as gene_exc
from genelocator.locate import BisectFinder


@pytest.fixture('module')
def build38finder():
    # Mostly, tests are slow because of the time to build this tree. It's immutable, so only do this 1x per run
    return get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True)


class TestGeneLocator:
    def test_finds_nearest_gene_when_none_overlap(self, build38finder):
        genes = build38finder.at('chr19', 1234)
        assert genes == [
            {'ensg': 'ENSG00000176695.8', 'chrom': '19', 'start': 107104, 'end': 117102, 'symbol': 'OR4F17'}
        ]

    def test_warns_on_unknown_chromosome(self, build38finder):
        with pytest.raises(gene_exc.BadCoordinateException, match="99"):
            build38finder.at('chr99', 1234)

    def test_finds_nearest_gene_in_region(self, build38finder):
        genes = build38finder.at('chr10', 112_950_250)
        assert genes == [
            {'ensg': 'ENSG00000148737.17', 'chrom': '10', 'start': 112_950_247, 'end': 113_167_678, 'symbol': 'TCF7L2'}
        ]

    def test_finds_multiple_genes_in_region(self, build38finder):
        genes = build38finder.at('10', 113_588_900)
        assert len(genes) == 2, 'This position overlaps two genes'
        assert genes == [{'chrom': '10',
                          'end': 113664127,
                          'ensg': 'ENSG00000197893.13',
                          'start': 113588716,
                          'symbol': 'NRAP'},
                         {'chrom': '10',
                          'end': 113589602,
                          'ensg': 'ENSG00000148702.15',
                          'start': 113550837,
                          'symbol': 'HABP2'}], 'Found two genes, and the nearest one was listed first'


class TestBisectFinder:
    """Validate that the bisect finder works as expected"""

    def test_scenarios(self):
        finder = BisectFinder([(23, 'foo'), (25, 'bar')])

        assert finder.get_item_before_or_at(22) is None, 'There is nothing before the beginning'
        assert finder.get_item_before_or_at(23) == (23, 'foo'), 'Found item matching target position'
        assert finder.get_item_before_or_at(24) == (23, 'foo'), 'Found nearest item before target position'
        assert finder.get_item_before_or_at(25) == (25, 'bar'), 'Found item matching target position'
        assert finder.get_item_before_or_at(26) == (25, 'bar'), 'Found nearest item before target position'

        assert finder.get_item_after_or_at(22) == (23, 'foo'), 'Found nearest item after target position'
        assert finder.get_item_after_or_at(23) == (23, 'foo'), 'Found nearest item at target position'
        assert finder.get_item_after_or_at(24) == (25, 'bar'), 'Found nearest item after target position'
        assert finder.get_item_after_or_at(25) == (25, 'bar'), 'Found nearest item at target position'
        assert finder.get_item_after_or_at(26) is None, 'There is nothing after the end'
