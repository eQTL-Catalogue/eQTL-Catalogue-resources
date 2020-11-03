# Tutorial on running all of the eQTL Catalogue workflow from start to finisih

Data processing for the eQTL Catalogue is based on the four main workflows:
* [eQTL-Catalogue/genimpute](https://github.com/eQTL-Catalogue/genimpute)
* [eQTL-Catalogue/rnaseq](https://github.com/eQTL-Catalogue/rnaseq)
* [eQTL-Catalogue/qcnorm](https://github.com/eQTL-Catalogue/qcnorm)
* [eQTL-Catalogue/qtlmap](https://github.com/eQTL-Catalogue/qtlmap)
* [eQTL-Catalogue/susie-workflow](https://github.com/eQTL-Catalogue/susie-workflow)

## Step 1: Genotype imputation with [eQTL-Catalogue/genimpute](https://github.com/eQTL-Catalogue/genimpute)

#### Input
Raw genotype data in PLINK binary format (.bed, .bim, .fam) using GRCh37 coordinates

#### Output
Imputed genotypes in VCF format lifted to GRCh38 coordinates.

#### Running the workflow
```bash
nextflow run main.nf -profile eqtl_catalogue -resume\
  --bfile <path_to_plink_file_prefix>\
  --harmonise_genotypes true\
  --output_name <output_file_prefix>\
  --outdir <path_to_output_dir>
```




