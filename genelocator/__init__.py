import gzip
import pickle

from genelocator.download import get_genes_iterator
from .locate import GeneLocator

from . import assets
from .const import BUILD_LOOKUP
from . import exception as gene_exc


def get_genelocator(build_or_path: str, *, gencode_version=32, coding_only=True, auto_fetch=False):
    if build_or_path in BUILD_LOOKUP:
        # We are looking up a special, known dataset cached on disk
        geneset = 'codinglike' if coding_only else 'all'  # TODO: Use enum here
        # If auto_fetch is specified, this function will block until the data has been returned
        # It would be really cool if we could use Peter's streaming dataset generation, but that saving wasn't
        # implemented. I've simplified the current code t reflect only the features we are currently using, with an eye
        # towards going that direction later.
        source_path = assets.locate_by_metadata(build_or_path, gencode_version, geneset, auto_fetch=auto_fetch)
    else:
        # The user has specified a path to a lookup file (premade tree in compressed pickle format). Use it!
        source_path = build_or_path

    with gzip.open(source_path, 'rb') as f:
        return pickle.load(f)
