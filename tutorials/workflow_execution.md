# Tutorial on running all eQTL Catalogue workflows from start to finisih

Data processing for the eQTL Catalogue is based on the four main workflows:
* [eQTL-Catalogue/genimpute](https://github.com/eQTL-Catalogue/genimpute)
* [eQTL-Catalogue/rnaseq](https://github.com/eQTL-Catalogue/rnaseq)
* [eQTL-Catalogue/qcnorm](https://github.com/eQTL-Catalogue/qcnorm)
* [eQTL-Catalogue/qtlmap](https://github.com/eQTL-Catalogue/qtlmap)

## Step 1: Genotype imputation with [eQTL-Catalogue/genimpute](https://github.com/eQTL-Catalogue/genimpute)

#### Input
Raw genotype data in PLINK binary format (.bed, .bim, .fam) using GRCh37 coordinates. If your data is in VCF format, then you need to first convert it to PLINK format with:

```bash
plink --vcf <path_to_vcf_file> --make-bed --out <plink_file_prefix>

```

Optionally, you can also immediately also check if there are some individual with many missing genotypes (see manual QC steps below). If there are then you should probably exclude those individuals, because their presence can cause the imputation workflow to fail.

```bash
plink --bfile <plink_file_prefix> --missing
```

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

#### Manual QC steps

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
3. Sample metadata file. See here for an example from the [GEUVADIS dataset](workflow_execution_files/GEUVADIS.tsv). Required columns:

#### Output
Normalised molecular trait (gene expression, exon expression, transcript usage, event usage) matrices in a format suitable for the qtlmap workflow.

#### Running the workflow
```bash
nextflow run main.nf\
  -profile tartu_hpc\
  --study_name <study_name>\
  --vcf_file <path_to_imputed_genotypes_from_the_genimpute_workflow.vcf.gz>\
  --exclude_population\
  --quant_results_path <rnaseq_workfow_output_folder_path>\
  --sample_meta_path <path_to_sample_metadata_file>\
  --skip_leafcutter_norm\
  --outdir <path_to_output_directory>\
  -process.queue amd\
  -resume
```

#### Manual QC steps

1. Check genotype concordance between the imputed genotypes and aligned RNA-seq reads. The summarised QTLtools mbv output is stored in the `QC/<study_name>_MBV_best_matches_matrix.tsv` file. Correct any obvious sample swaps between RNA-seq and genotype files. Exclude RNA-seq samples that do not match any individuals in the genotype data by setting the `genotype_qc_passed` field to FALSE in the sample metadata file. Also exclude RNA-seq samples that match more than one individual in the genotype data as this could be a sign of cross-contamination between RNA samples (NOTE: make sure that each individual is present only once in the genotype data and your data does not contain any monoxygotic twins.). 
2. Study the 'pop_assing/relatedness_matrix.tsv' file to ensure that there are no related individuals in the genotype data. Relatedness values > 0.2 are an obvious concern.
3. Check the gene expression PCA and MDS plots in the `QC/<study_name>_QC_report.html` file. Exclude any obvious outliers by setting `rna_qc_passed` field to FALSE in the sample metadata file. 
4. Check of the expression of sex-specific genes is consistent with the annotated sex of the samples. Fix missing or mis-annotated sex in the sample metadata file. Note that high simultaneous expression of both XIST (female-specifc) and Y chromosome genes (male-specifc) can be a good indiciation of RNA cross-contamination between two samples. This is often concordant with the results seen in the RNA-seq analysis.


## Step 4: QTL analysis and fine mapping with [eQTL-Catalogue/qtlmap](https://github.com/eQTL-Catalogue/qtlmap)

#### Input
1. Normalised moleocular trait files from the qcnorm workflow. **NOTE:** The studyFile from qcnorm `(qcnorm_output_directory>/<study_name>/<study_name>_qtlmap_inputs.tsv)` contains relative paths. You should either copy the qcnorm output directory to the qtlmap directory or create a symlink with the same name. 

#### Output

#### Running the workflow

```bash
nextflow run main.nf -profile eqtl_catalogue\
  --studyFile <qcnorm_output_directory>/<study_name>/<study_name>_qtlmap_inputs.tsv\
  --vcf_has_R2_field FALSE\
  --varid_rsid_map_file testdata/varid_rsid_map.tsv.gz\
  --n_batches 200
```


