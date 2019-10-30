"""
A command-line script that identifies nearest gene based on genomic coordinates

Sample usages:
    gene-locator GRCh37 gencode32 codinglike-genes chr19 234523
    gene-locator hg38 gencode31 all-genes 18 234523
"""

import argparse

from genelocator import get_genelocator


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
                        choices=["GRCh37", "GRCh38", "hg19", "hg38"],
                        help="The genome build")
    parser.add_argument("version", default=32, type=_validate_gencode,
                        help="The GENCODE database version to use")
    parser.add_argument("gene_type", choices=["codinglike-genes", "all-genes"],
                        # TODO: Do we expect more than 2 options? If not, a boolean flag might be more appropriate here
                        # TODO: Rename values
                        help='The type of gene: "codinglike-genes" (which means protein_coding_gene + IG_*_gene + TR_*_gene) or "all-genes"')  # noqa
    parser.add_argument("chromosome")
    parser.add_argument("position", type=int,
                        help="The variant position (coordinates should refer to the selected genome build")
    return parser.parse_args()


def main():
    args = parse_args()
    build = args.build

    # TODO don't coerce to int?
    grch_build = int(build[-2:])  # Just numeric part (this is a bit hacky)
    if grch_build == 19:
        grch_build = 37

    gencode_version = args.version
    coding_only = (args.gene_type == 'codinglike')
    chrom = args.chromosome.replace('chr', '')
    genelocator = get_genelocator(grch_build, gencode_version, coding_only)
    for gene in genelocator.at(chrom, args.position):
        print('{chrom}\t{start}\t{end}\t{ensg}\t{symbol}'.format(**gene))


if __name__ == '__main__':
    main()
