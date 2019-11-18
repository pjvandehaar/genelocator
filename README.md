## Gene Locator
This library allows the user to find which genes overlap a genomic location, given a chromosome and position.
If no genes overlap the given position, then the closest gene is returned instead.

### Usage
This library can be used via command line script, or called from within Python. Only Python versions >=3.5
are supported.

```sh
$ pip3 install genelocator
$ gene-locator GRCh37 chr19 234523 --common-genetypes --version gencode32
# => 19	281040	291403	ENSG00000141934.10_5	PLPP2
```

```python3
from genelocator import get_genelocator
# By default, it will only perform the lookup if cached data is available.
#  A new lookup can be automatically generated for a different build/ gene list, by specifying auto_fetch=True
gl = get_genelocator('GRCh38', gencode_version=31, common_genetypes_only=True, auto_fetch=True)
gene = gl.at('chr19', 101000)
# => [{'chrom': '19', 'start': 107104, 'end': 117102, 'ensg': 'ENSG00000176695.8', 'symbol': 'OR4F17'}]
```

The python package comes bundled with data from GENCODE version 32, for builds GRCh37 and GRCh38.


### Rules
It works as follows:

1. If the genomic location (chromosome and position) falls within at least one gene, return information about each overlapped gene.
2. If the location does not fall within any genes, return information for the gene whose start (TSS) or end (TES) is closest.
3. If the requested chromosome has no data, throw an error.


### Development
To install dependencies and run in development mode:

`pip install -e .`

Linting rules: `flake8`

Unit tests: `pytest tests/`
