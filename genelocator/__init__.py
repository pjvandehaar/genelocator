
from .download import get_genes
from .locate import GeneLocator


def get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    return GeneLocator(get_genes(grch_build, gencode_version, only_codinglike_genetypes))
