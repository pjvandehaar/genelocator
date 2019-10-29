
def test_all():
    from genelocator import get_genelocator
    gl = get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True)
    gene = gl.at('chr19', 1234)
    assert gene == [('ENSG00000176695.8', '19', 107104, 117102, 'OR4F17')]
