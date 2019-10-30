"""
This script downloads Gencode data and converts

This module exists purely to provide get_genes()
"""
import gzip
import io
import json
import os
import re
import urllib.request


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


# FIXME: Add an autofetch option: some use cases (like web apps) don't want to silently block the first time someone
#  requests an unknown build
def get_genes(*,
              grch_build: int = 38,
              gencode_version: int = 32,
              only_codinglike_genetypes: bool = True):
    """
    Get a list of genes (represented by dicts). The CODINGLIKE_GENETYPES in this module were chosen manually.

    This function well attempt to generate the lookup locally,

    """
    filepath = os.path.join(os.path.dirname(__file__),
                            'data',
                            _get_filename(grch_build=grch_build, gencode_version=gencode_version,
                                          only_codinglike_genetypes=only_codinglike_genetypes))
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
                if only_codinglike_genetypes and genetype not in CODINGLIKE_GENETYPES:
                    continue
                yield {'chrom': chrom, 'start': start, 'end': end, 'ensg': ensg, 'symbol': symbol}


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


def _get_filename(*, grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    return 'gencode-grch{}-gencode{}-{}.json.gz'.format(grch_build, gencode_version,
                                                        'codinglike' if only_codinglike_genetypes else 'all')


def _save_file(grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    """cache get_genes() results into a file in this directory for quick use later"""
    genes = list(get_genes(
        grch_build=grch_build,
        gencode_version=gencode_version,
        only_codinglike_genetypes=only_codinglike_genetypes
    ))
    compressed = gzip.compress(json.dumps(genes, separators=(',', ':')).encode('utf8'))
    filename = _get_filename(grch_build=grch_build,
                             gencode_version=gencode_version,
                             only_codinglike_genetypes=only_codinglike_genetypes)
    with open(filename, 'wb') as f:
        f.write(compressed)


def _save_all_files(gencode_version=32):
    """cache all data for a single (ideally latest) gencode version"""
    for grch_build in [37, 38]:
        _save_file(grch_build, gencode_version)
        _save_file(grch_build, gencode_version, False)
