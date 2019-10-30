## Gene Locator
This library allows the user to annotate a list of genetic variants given chromosome/position as coordinates.

### Usage
This library can be used via command line script, or called from within Python. Currently only Python versions >=3.5 
are supported.

```sh
$ pip3 install genelocator
$ gene-locator GRCh37 chr19 234523 --coding-only --version gencode32
# => 19	281040	291403	ENSG00000141934.10_5	PLPP2
```

```python3
from genelocator import get_genelocator
# By default, it will only perform the lookup if cached data is available.  
#  A new lookup can be automatically generated for a different build/ gene list, by specifying auto_fetch=True
gl = get_genelocator('GRCh38', gencode_version=31, coding_only=True, auto_fetch=True)
gene = gl.at('chr19', 101000)
# => [{'chrom': '19', 'start': 107104, 'end': 117102, 'ensg': 'ENSG00000176695.8', 'symbol': 'OR4F17'}]
```

The python package comes bundled with data from GENCODE version 32, for builds GRCh37 and GRCh38. 


### Rules
It works as follows:
1. If a SNP falls within at least one gene, return a list of gene information for each gene
    1a. If a SNP falls within multiple genes, return a list of information about all overlapping genes, sorted by
     nearest first
2. If a SNP does not fall within any genes, return information for the gene whose start or end is closest 
to the specified coordinates.
3. If the requested chromosome has no data, throw an error.


### Development
To install dependencies and run in development mode:

`pip install -e .`

Linting rules: `flake8`

Unit tests: `pytest tests/` 
