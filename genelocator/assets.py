"""
Contains helpers that can be used to locate (or generate) asset data files

For now, assets will be cached in the `data` folder. Eventually we can extract this to a package.
"""
import logging
import os

from .const import BUILD_LOOKUP, KNOWN_GENESETS
from . import download
from . import exception as gene_exc

logger = logging.getLogger(__name__)


def _get_cache_filepath(build: str, version: int, geneset: str) -> str:
    """Get the path to a cached, pre-built copy of the asset"""
    build_numeric = BUILD_LOOKUP[build]
    filename = 'genes-grch{}-gencode{}-{}.pickle.gz'.format(build_numeric, version, geneset)
    return os.path.join(os.path.dirname(__file__), 'data', filename)


def _get_generated_filepath(build: str, version: int, geneset: str) -> str:
    """
    Trigger generating/fetching the required dataset (if possible). This will block until completed, so it may be slow.

    Generate a dataset (if possible). Returns the filename of the generated asset, or raises an AssetFetchError
    """
    out_path = _get_cache_filepath(build, version, geneset)
    coding_only = (geneset == 'codinglike')
    download.make_gene_locator(build, out_path,
                               gencode_version=version, coding_only=coding_only)
    return out_path


def locate_by_metadata(build: str, version: int, geneset: str, *, auto_fetch=False) -> str:
    """Locate an asset file that satisfies all specified parameters"""
    if geneset not in KNOWN_GENESETS:
        # Builds or gencode versions might be a new file (in which case, checking from a server is ok)
        # But how genes are selected is a set of rules defined in code, so these labels are rigidly defined
        raise gene_exc.UnsupportedDatasetException

    # Best option: find a cached copy of the asset on disk
    target_filename = _get_cache_filepath(build, version, geneset)
    if os.path.isfile(target_filename):
        return target_filename
    elif not auto_fetch:
        raise gene_exc.NoCachedDataException('Failed to locate requested dataset: {}'.format(target_filename))

    # There is no cached copy of the dataset, and the user has chosen to auto-fetch a copy
    # This function will download (or generate) the asset, and return a path when the process has completed
    logger.info("No cached asset found; attempting to download")
    return _get_generated_filepath(build, version, geneset)
