
import sys
from . import get_genelocator


def main():
    if len(sys.argv) != 6:
        print("This command prints information about the closest gene (or genes if they all overlap the location).")
        print("")
        print("Usage:")
        print("    genelocator <genome_build> <gencode_version> <genetypes> <chromosome> <position>")
        print("Examples:")
        print("    genelocator GRCh37 gencode32 codinglike-genes chr19 234523")
        print("    genelocator hg38 gencode31 all-genes 18 234523")
        sys.exit(1)
    if sys.argv[1].lower() in ('grch37', 'grch38', 'hg19', 'hg38'):
        grch_build = int(sys.argv[1][-2:])
        if grch_build == 19: grch_build = 37
    else:
        print('Genome build must be "GRCh37" or "GRCh38", not {!r}'.format(sys.argv[1]))
        sys.exit(1)
    if sys.argv[2].lower().startswith('gencode') and sys.argv[2][7:].isdigit():
        gencode_version = int(sys.argv[2][7:])
        if gencode_version <= 22:
            print("Gencode versions 22 and earlier are not supported")
            sys.exit(1)
    if sys.argv[3].lower().startswith(('codinglike', 'all')):
        only_codinglike_genetypes = sys.argv[3].lower().startswith('codinglike')
    else:
        print('Genetypes must be "codinglike-genes" (which means protein_coding_gene + IG_*_gene + TR_*_gene) or "all-genes", not {!r}'.format(sys.argv[3]))
        sys.exit(1)
    chrom = sys.argv[4] if not sys.argv[4].startswith('chr') else sys.argv[4][3:]
    pos = int(sys.argv[5])
    genelocator = get_genelocator(grch_build, gencode_version, only_codinglike_genetypes)
    for gene in genelocator.at(chrom, pos):
        print('{chrom}\t{start}\t{end}\t{ensg}\t{symbol}'.format(**gene))
