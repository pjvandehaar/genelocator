## Gene Locator
This library allows the user to annotate a list of genetic variants given chromosome/position as coordinates.

### Usage
This library can be used via command line script, or called from within Python. Currently only Python versions >=3.5 
are supported.

```sh
$ pip3 install genelocator
$ gene-locator GRCh37 gencode32 codinglike-genes chr19 290123
# => 19	281040	291403	ENSG00000141934.10_5	PLPP2
```

```python3
from genelocator import get_genelocator
gl = get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True)
gene = gl.at('chr19', 101000)
# => [{'chrom': '19', 'start': 107104, 'end': 117102, 'ensg': 'ENSG00000176695.8', 'symbol': 'OR4F17'}]
```


### Rules
It works as follows:
1. If a SNP falls within at least one gene, return a list of gene information for each gene
    1a. If a SNP falls within multiple genes, return _____ (how sorted?)
2. If a SNP does not fall within any genes, return information for the gene whose transcription start site is closest 
to the specified coordinates.
3. If the requested chromosome has no data, ____.


### Development
To install dependencies and run in development mode:

`pip install -e .`

Linting rules: `flake8`

Unit tests: `pytest tests/` 
