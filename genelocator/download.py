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


# These "common" genetypes are usually the most useful
# To see all genetypes, run `Counter(g['genetype'] for g in _get_unfiltered_genes_iterator(38)).most_common()`
# These genetypes are copied from <https://github.com/hyunminkang/cramore/blob/6d85ed9/cmd_plp_make_dge_matrix.cpp#L137>.
COMMON_GENETYPES = {
    'protein_coding',
    'IG_C_gene',
    'IG_D_gene',
    'IG_J_gene',
    'IG_V_gene',
    'TR_C_gene',
    'TR_D_gene',
    'TR_J_gene',
    'TR_V_gene'
    'lincRNA',
    'Mt_tRNA',
    'Mt_rRNA',
    'antisense',
}


def _download_gencode_gtf_gz_bytes(grch_build_number: int, gencode_version: int) -> bytes:
    """Download the file, hold it in memory, and return it as bytes"""
    if grch_build_number == 37:
        template = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/GRCh37_mapping/gencode.v{gencode_version}lift37.annotation.gtf.gz'
    elif grch_build_number == 38:
        template = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/gencode.v{gencode_version}.basic.annotation.gtf.gz'
    else:
        raise Exception('cannot handle GRCh build {!r}'.format(grch_build_number))
    with urllib.request.urlopen(url=template.format(gencode_version=gencode_version)) as f:
        x = f.read()
        content_length = int(f.headers.get('Content-length'))
        if len(x) != content_length:
            raise Exception('Expected {!r} bytes but received {!r} bytes beginning {!r}.'.format(
                content_length, len(x), x[:100]))
        return x


def _get_genelist_filename(grch_build_number: int, *, gencode_version: int = 32, common_genetypes_only: bool = True) -> str:
    return 'gencode-{}-gencode{}-{}.json.gz'.format(grch_build_number, gencode_version,
                                                    'common_genetypes_only' if common_genetypes_only else 'all')


def _save_locator(locator: GeneLocator, out_path):
    """
    Once the interval tree has been created, it is convenient to save it for future use.
    """
    with gzip.open(out_path, 'wb') as f:
        # protocol=4 is faster and supported in python3.4+
        pickle.dump(locator, f, protocol=4)


def _get_unfiltered_genes_iterator(grch_build_number: int, *, gencode_version: int = 32) -> ty.Iterator[dict]:
    compressed = _download_gencode_gtf_gz_bytes(grch_build_number, gencode_version)
    with gzip.open(io.BytesIO(compressed), 'rt') as f:  # avoid holding the full text uncompressed in RAM at once
        for line in f:
            if line.startswith('#'):
                continue
            chrom, source, feature_type, start, end, _, strand, CDS_phase, info = line.split('\t')
            if feature_type != 'gene':
                continue
            start = int(start)
            end = int(end)
            if start >= end:
                raise Exception('start >= end for line {!r}'.format(line))
            ensg = _extract_first_group(r'gene_id "(ENSGR?[0-9._A-Z]+?)"', info)  # Sometimes we want `ensg.split('.')[0]` but not here.
            symbol = _extract_first_group(r'gene_name "(.+?)"', info)
            genetype = _extract_first_group(r'gene_type "(.+?)"', info)
            yield {'chrom': chrom, 'start': start, 'end': end, 'ensg': ensg, 'symbol': symbol, 'genetype': genetype, 'line': line}


def get_genes_iterator(grch_build: str, *, gencode_version: int = 32, common_genetypes_only: bool = True) -> ty.Iterator[dict]:
    """
    Get a list of genes (represented by dicts). The CODINGLIKE_GENETYPES in this module were chosen manually.
    This creates an intermediate datafile that is then used to build the interval tree, so usually this function is
        only used during initial dataset generation
        (not during day-to-day analysis)
    """
    try:
        grch_build_number = const.BUILD_LOOKUP[grch_build]
    except KeyError:
        raise gene_exc.AssetFetchError('Cannot retrieve assets for build {}'.format(grch_build))

    filepath = os.path.join(os.path.dirname(__file__),
                            'data',
                            _get_genelist_filename(grch_build_number, gencode_version=gencode_version,
                                                   common_genetypes_only=common_genetypes_only))
    if os.path.exists(filepath):
        with gzip.open(filepath, 'rt') as f:
            yield from json.load(f)
    else:
        for gene in _get_unfiltered_genes_iterator(grch_build_number, gencode_version=gencode_version):
            if not gene['chrom'].startswith('chr'):
                if gene['chrom'].startswith('GL'):
                    # We don't know what to do with this type of entry, even though it's a valid identifier
                    continue
                else:
                    raise Exception('Unknown chromosome {!r} on line {!r}'.format(gene['chrom'], gene['line']))
            if common_genetypes_only and gene['genetype'] not in COMMON_GENETYPES:
                continue
            gene.pop('genetype')
            gene.pop('line')
            yield gene


def _extract_first_group(pattern: str, string: str) -> str:
    match = re.search(pattern, string)
    if match is None:
        raise Exception("re.search({!r}, {!r}) didn't match!")
    return match.group(1)


def make_gene_locator(grch_build: str, out_path, *, gencode_version: int = 32, common_genetypes_only: bool = True) -> GeneLocator:
    # TODO: We should save this genes list for performance reasons if possible
    genes = get_genes_iterator(grch_build, gencode_version=gencode_version, common_genetypes_only=common_genetypes_only)
    locator = GeneLocator(genes)
    _save_locator(locator, out_path)
    return locator
