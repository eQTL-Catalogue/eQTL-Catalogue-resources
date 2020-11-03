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
  --output_name <output_prefix>\
  --outdir <path_to_output_dir>
```

#### Manual quality control steps

- Check the <output_prefix>.imiss file for samples with large proportion of missing genotypes (e.g. > 5%). These samples are likely to have poor genotyping quality and should probably the excluded from the analysis before continuing. Remove these inviduals from the original plink file and re-run the genimpute workflow.

## Step 2: RNA-seq quantification with [eQTL-Catalogue/rnaseq](https://github.com/eQTL-Catalogue/rnaseq)

## Step 3: Gene expression and genotype data normalisation and QC with [eQTL-Catalogue/qcnorm](https://github.com/eQTL-Catalogue/qcnorm)

## Step 4: QTL mapping with [eQTL-Catalogue/qtlmap](https://github.com/eQTL-Catalogue/qtlmap)

## Step 5: QTL fine mapping with [eQTL-Catalogue/susie-workflow](https://github.com/eQTL-Catalogue/susie-workflow)

# Modular approach



