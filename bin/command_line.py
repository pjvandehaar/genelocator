"""
A command-line script that identifies nearest gene based on genomic coordinates

Sample usages:
    gene-locator GRCh37 chr19 234523 --coding-only --version gencode32
    gene-locator hg38 18 234523 --version gencode31
"""

import argparse
import logging
import sys

from genelocator import get_genelocator
from genelocator.const import BUILD_LOOKUP
from genelocator import exception as gene_exc


logger = logging.getLogger(__name__)


def _validate_gencode(value) -> int:
    """Due to file format changes, we don't support older data versions"""
    # Strip name prefix if present
    if value.startswith('gencode') and value[7:].isdigit():
        value = int(value[7:])

    value = int(value)
    if value <= 22:
        raise argparse.ArgumentTypeError("Gencode versions 22 and earlier are not supported")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="This command prints information about the closest gene (or genes, if multiple overlap the location).")  # noqa
    parser.add_argument("build",
                        choices=BUILD_LOOKUP.keys(),
                        help="The genome build (must be specified)")
    parser.add_argument("chromosome")
    parser.add_argument("position", type=int,
                        help="The variant position (coordinates should refer to the selected genome build")
    parser.add_argument("--version", default=32, type=_validate_gencode,
                        help="The GENCODE database version to use")
    parser.add_argument("--coding-only", dest="coding_only", action='store_true',
                        help='If specified, restrict search to "coding-like" genes (which means protein_coding_gene + IG_*_gene + TR_*_gene)')  # noqa
    parser.add_argument('--auto-fetch', dest='auto_fetch', action='store_true',
                        help="If specified, will automatically try to download the required data")
    return parser.parse_args()


def main():
    # TODO (future) Creating the tree is the slow step. A potential performance enhancement would be to allow searching
    #   for multiple positions from a single run of the command line (so we get many searches but pay startup penalty once)
    args = parse_args()
    gencode_version = args.version
    try:
        genelocator = get_genelocator(args.build,
                                      gencode_version=gencode_version,
                                      coding_only=args.coding_only,
                                      auto_fetch=args.auto_fetch)
    except gene_exc.UnsupportedDatasetException:
        logger.error('No source found for the requested dataset; exiting')
        sys.exit(1)

    for gene in genelocator.at(args.chromosome, args.position):
        print('{chrom}\t{start}\t{end}\t{ensg}\t{symbol}'.format(**gene))


if __name__ == '__main__':
    main()
