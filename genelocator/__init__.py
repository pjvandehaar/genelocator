import gzip
import pickle

from genelocator.download import get_genes_iterator  # noqa: F401
from .locate import GeneLocator  # noqa: F401

from . import assets
from .const import BUILD_LOOKUP
from . import exception as gene_exc  # noqa: F401


def get_genelocator(build_or_path: str, *, gencode_version: int = 32, common_genetypes_only=True, coding_only=None, auto_fetch=False) -> GeneLocator:

    # for backward-compatibility with old argument name
    if coding_only is not None:
        common_genetypes_only = coding_only

    if build_or_path in BUILD_LOOKUP:
        # We are looking up a special, known dataset cached on disk
        geneset = 'common_genetypes' if common_genetypes_only else 'all'  # TODO: Use enum here
        # If auto_fetch is specified, this function will block until the data has been returned
        source_path = assets.locate_by_metadata(build_or_path, gencode_version, geneset, auto_fetch=auto_fetch)
    else:
        # The user has specified a path to a lookup file (premade tree in compressed pickle format). Use it!
        source_path = build_or_path

    with gzip.open(source_path, 'rb') as f:
        return pickle.load(f)
