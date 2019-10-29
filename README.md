
```sh
pip3 install genelocator
genelocator GRCh37 gencode32 codinglike-genes chr19 290123
# => 19	281040	291403	ENSG00000141934.10_5	PLPP2
```

```python3
from genelocator import get_genelocator
gl = get_genelocator(grch_build=38, gencode_version=32, only_codinglike_genetypes=True)
gene = gl.at('chr19', 101000)
# => [{'chrom': '19', 'start': 107104, 'end': 117102, 'ensg': 'ENSG00000176695.8', 'symbol': 'OR4F17'}]
```
