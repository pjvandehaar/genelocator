

def get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    import pickle, os
    from .download import get_genes
    from .locate import GeneLocator
    filename = _get_filename(grch_build, gencode_version, only_codinglike_genetypes)
    filepath = os.path.join(os.path.dirname(__file__), filename)
    if os.path.exists(filepath):
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    else:
        return GeneLocator(get_genes(grch_build, gencode_version, only_codinglike_genetypes))


def _get_filename(grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    return 'gencode-grch{}-gencode{}-{}.pickle'.format(grch_build, gencode_version, 'codinglike' if only_codinglike_genetypes else 'all')

def _save_file(grch_build=38, gencode_version=32, only_codinglike_genetypes=True):
    '''cache get_genelocator() result into a file in this directory for quick use later'''
    import pickle
    gl = get_genelocator(grch_build, gencode_version, only_codinglike_genetypes)
    filename = _get_filename(grch_build, gencode_version, only_codinglike_genetypes)
    with open(filename, 'wb') as f:
        # protocol=4 is faster and supported in python3.4+
        pickle.dump(gl, f, protocol=4)
