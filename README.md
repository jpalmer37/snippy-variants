[![Tests](https://github.com/BCCDC-PHL/snippy-variants/actions/workflows/tests.yml/badge.svg)](https://github.com/BCCDC-PHL/snippy-variants/actions/workflows/tests.yml)

# snippy-variants

Perform read mapping and variant calling using [snippy](https://github.com/tseemann/snippy).

## Usage

```
nextflow run BCCDC-PHL/snippy-variants \
  --ref </path/to/ref.fa> \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

