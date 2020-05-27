# Explanation of the columns in the eQTL summary statistics files

* **variant** - The variant ID (chromosome_position_ref_alt) e.g. 19_226776_C_T. Based on GRCh38 coordinates and reference genome. The chromosome, position, ref and alt values should exactly match same fields in the summary statistics file. 
* **r2** - Optional imputation quality score from the imputation software, can be replaced with NA if not available.
* **pvalue** - Nominal p-value of association between the variant and the molecular trait.
* **molecular_trait_object_id** - For phenotypes with multiple correlated alternatives (multiple alternative transcripts or exons within a gene, multple alternative promoters in txrevise, multiple alternative intons in Leafcutter), this defines the level at which the phenotypes were aggregated. Permutation p-values are calculated accross this set of alternatives.  
* **molecular_trait_id** - ID of the molecular trait used for QTL mapping. Depending on the quantification method used, this can be either a gene id, exon id, transcript id or a txrevise promoter, splicing or 3'end event id. Examples: ENST00000356937, ENSG00000008128.  
* **maf** - Minor allele frequency within a QTL mapping context (e.g. cell type or tissues within a study).
* **gene_id** - Ensembl gene ID of the molecular trait. 
* **median_tpm** - Median transcripts per million (TPM) expression value of the gene. Can be replaced with NA if not availble (e.g. in microarray  studies).
* **beta** - Regression coefficient from the linear model.
* **se** - Standard error of the beta.
* **an** - Total number of alleles. For autosomal variants, this is usually two times the sample size. Conversly, for autosomal variants, sample size is equal to an/2.
* **ac** - Count of the alternative allele. 
* **chromosome** - GRCh38 chromosome name of the variant (e.g. 1,2,3 ...,X).
* **position** - GRCh38 position of the variant.
* **ref** - GRCh38 reference allele.
* **alt** - GRCh38 alternative allele (also the effect allele).
* **type** - Type of the genetic variant; SNP, INDEL or OTHER.
* **rsid** - The dbSNP v151 rsid of the variant. If the same variant has multiple rsids then these should be split over multiple rows so that all of the other values are duplicates.
