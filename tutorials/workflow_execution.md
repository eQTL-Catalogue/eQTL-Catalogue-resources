# Tutorial on running all eQTL Catalogue workflows from start to finisih

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

#### Input
1. Raw RNA-seq data in fastq format. 
  - Is it paired-end or single-end? (Do you have one or two fastq files per sample?)
  - Is it stranded or unstranded?
  - Paths to the fastq files can be passed with the readPathsFile parameter. See example files for [paired-end](https://github.com/eQTL-Catalogue/rnaseq/blob/master/data/readPathsFile_macrophages_PE.tsv) and [single-end](https://github.com/eQTL-Catalogue/rnaseq/blob/master/data/readPathsFile_macrophages_SE.tsv) data.
2. Imputed genotypes in VCF format (from the genimpute workflow). These are required to check genotype condordance between the VCF files and the RNA-seq data using the [MBV](https://doi.org/10.1093/bioinformatics/btx074) method.

#### Output
Raw gene expression, exon expression, transcript expression and event expression matrices in a format suitable for the eQTL-Catalogue/qcnorm workflow.

#### Running the workflow
```bash
nextflow run main.nf\
 -profile eqtl_catalogue\
 --readPathsFile <path_to_readPathsFile.tsv>\
 --reverse_stranded\
 --skip_qc\
 --skip_multiqc\
 --skip_stringtie\
 --saveReference\
 --run_tx_exp_quant\
 --run_txrevise\
 --run_splicing_exp_quant\
 --run_exon_quant\
 --run_mbv\
 --mbv_vcf <path_to_imputed_genotypes_from_the_genimpute_workflow.vcf.gz>\
 -process.queue main\
 -resume
```

#### Other useful options
Use the `-executor.queueSize` option to limit the number alignment jobs running in parallel to avoid too much load on the disks.


## Step 3: Gene expression and genotype data normalisation and QC with [eQTL-Catalogue/qcnorm](https://github.com/eQTL-Catalogue/qcnorm)

#### Input
1. Imputed genotypes from the genimpute workflow
2. RNA-seq quantification results from the rnaseq workflow.
3. Sample metadata file. See here for an example from the GEUVADIS dataset. Required columns:

#### Output
Normalised molecular trait (gene expression, exon expression, transcript usage, event usage) matrices in a format suitable for the qtlmap workflow.

```bash
nextflow run main.nf\
  -profile tartu_hpc\
  --study_name <name_of_the_study>\
  --vcf_file <path_to_imputed_genotypes_from_the_genimpute_workflow.vcf.gz>\
  --exclude_population\
  --quant_results_path <rnaseq_workfow_output_folder_path>\
  --sample_meta_path <path_to_sample_metadata_file>\
  --skip_leafcutter_norm\
  --outdir <path_to_output_directory>\
  -process.queue amd\
  -resume
```

## Step 4: QTL mapping with [eQTL-Catalogue/qtlmap](https://github.com/eQTL-Catalogue/qtlmap)

## Step 5: QTL fine mapping with [eQTL-Catalogue/susie-workflow](https://github.com/eQTL-Catalogue/susie-workflow)

# Modular approach



