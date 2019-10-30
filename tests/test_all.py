from genelocator import get_genelocator
from genelocator.locate import BisectFinder


def test_lookup():
    gl = get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True)
    genes = gl.at('chr19', 1234)
    assert genes == [
        {'ensg': 'ENSG00000176695.8', 'chrom': '19', 'start': 107104, 'end': 117102, 'symbol': 'OR4F17'}
    ]


class TestBisectFinder:
    def test_scenarios(self):
        """These tests were predefined. TODO: Add better explanation of what each scenario is testing"""
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_before_or_at(22) is None
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_before_or_at(23) == (23, 'foo')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_before_or_at(24) == (23, 'foo')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_before_or_at(25) == (25, 'bar')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_before_or_at(26) == (25, 'bar')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_after_or_at(22) == (23, 'foo')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_after_or_at(23) == (23, 'foo')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_after_or_at(24) == (25, 'bar')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_after_or_at(25) == (25, 'bar')
        assert BisectFinder([(23, 'foo'), (25, 'bar')]).get_item_after_or_at(26) is None
