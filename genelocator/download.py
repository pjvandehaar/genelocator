"""
This script downloads Gencode data and converts it to a locator object

This represents the "build step" to create a locator object from raw data
"""
import gzip
import io
import json
import os
import pickle
import re
import typing as ty
import urllib.request

from . import const
from . import exception as gene_exc
from .locate import GeneLocator


# These "coding-like" genetypes are usually the most useful
CODINGLIKE_GENETYPES = {
    'protein_coding',
    'IG_C_gene',
    'IG_D_gene',
    'IG_J_gene',
    'IG_V_gene',
    'TR_C_gene',
    'TR_D_gene',
    'TR_J_gene',
    'TR_V_gene'
}


def _download_gencode_gtf_gz_bytes(grch_build, gencode_version):
    """Download the file, hold it in memory, and return it as bytes"""
    if grch_build == 37:
        template = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/GRCh37_mapping/gencode.v{gencode_version}lift37.annotation.gtf.gz'  # noqa
    elif grch_build == 38:
        template = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/gencode.v{gencode_version}.basic.annotation.gtf.gz'  # noqa
    else:
        raise Exception('cannot handle GRCh build {!r}'.format(grch_build))
    with urllib.request.urlopen(url=template.format(gencode_version=gencode_version)) as f:
        x = f.read()
        content_length = int(f.headers.get('Content-length'))
        if len(x) != content_length:
            raise Exception('Expected {!r} bytes but received {!r} bytes beginning {!r}.'.format(
                content_length, len(x), x[:100]))
        return x


def _get_genelist_filename(grch_build: str, *, gencode_version: int = 32, coding_only: bool = True):
    return 'gencode-{}-gencode{}-{}.json.gz'.format(grch_build, gencode_version,
                                                    'codinglike' if coding_only else 'all')


def _save_locator(locator: GeneLocator, out_path: str):
    """
    Once the interval tree has been created, it is convenient to save it for future use.
    """
    with gzip.open(out_path, 'wb') as f:
        # protocol=4 is faster and supported in python3.4+
        pickle.dump(locator, f, protocol=4)


def get_genes_iterator(grch_build: str, *, gencode_version: int = 32, coding_only: bool = True) -> ty.Iterator:
    """
    Get a list of genes (represented by dicts). The CODINGLIKE_GENETYPES in this module were chosen manually.
    This creates an intermediate datafile that is then used to build the interval tree, so usually this function is
        only used during initial dataset generation
        (not during day-to-day analysis)
    """
    try:
        grch_build = const.BUILD_LOOKUP[grch_build]
    except KeyError:
        raise gene_exc.AssetFetchError('Cannot retrieve assets for build {}'.format(grch_build))

    filepath = os.path.join(os.path.dirname(__file__),
                            'data',
                            _get_genelist_filename(grch_build, gencode_version=gencode_version,
                                                   coding_only=coding_only))
    if os.path.exists(filepath):
        with gzip.open(filepath, 'rt') as f:
            yield from json.load(f)
    else:
        compressed = _download_gencode_gtf_gz_bytes(grch_build, gencode_version)
        with gzip.open(io.BytesIO(compressed), 'rt') as f:  # avoid holding the full text uncompressed in RAM at once
            for line in f:
                if line.startswith('#'):
                    continue
                chrom, source, feature_type, start, end, _, strand, CDS_phase, info = line.split('\t')
                if feature_type != 'gene':
                    continue
                start = int(start)
                end = int(end)
                if not chrom.startswith('chr'):
                    if chrom.startswith('GL'):
                        # We don't know what to do with this type of entry, even though it's a valid identifier
                        continue
                    else:
                        raise Exception('Unknown chromosome {!r} on line {!r}'.format(chrom, line))
                assert start < end, line
                ensg = re.search(r'gene_id "(ENSGR?[0-9._A-Z]+?)"', info).group(
                    1)  # Sometimes we want `ensg.split('.')[0]` but not here.
                symbol = re.search(r'gene_name "(.+?)"', info).group(1)
                genetype = re.search(r'gene_type "(.+?)"', info).group(1)
                if coding_only and genetype not in CODINGLIKE_GENETYPES:
                    continue
                yield {'chrom': chrom, 'start': start, 'end': end, 'ensg': ensg, 'symbol': symbol}


def make_gene_locator(grch_build: str, out_path, *, gencode_version: int = 32, coding_only: bool = True):
    # TODO: We should save this genes list for performance reasons if possible
    genes = get_genes_iterator(grch_build, gencode_version=gencode_version, coding_only=coding_only)
    locator = GeneLocator(genes)
    _save_locator(locator, out_path)
    return locator

# def _save_genes_list(grch_build: str, *, gencode_version=32, coding_only=True):
#     """cache get_genes_iterator() results into a file in this directory for quick use later"""
#     genes = list(get_genes_iterator(grch_build=grch_build, gencode_version=gencode_version, coding_only=coding_only))
#     compressed = gzip.compress(json.dumps(genes, separators=(',', ':')).encode('utf8'))
#     filename = _get_genelist_filename(grch_build, gencode_version=gencode_version, coding_only=coding_only)
#     with open(filename, 'wb') as f:
#         f.write(compressed)


if __name__ == '__main__':
    pass
    # gencode_version = 32
    # for grch_build in [37, 38]:
    #     _save_genes_list(grch_build, gencode_version=gencode_version)
    #     _save_genes_list(grch_build, gencode_version=gencode_version, coding_only=False)
