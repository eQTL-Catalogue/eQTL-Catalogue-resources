# Files available from the FTP server

* **QTD\*.all.parquet** - Full nominal cis-QTL summary staistics for all tested molecular traits (available only for ge, microarray and aptamer quantification methods).
* **QTD\*.cc.parquet** - Nominal cis-QTL summary staistics for all molecular traits that permutation-based FDR < 1% and at least one fine mapped credible set detected by SuSiE. This is a filtered version of the **QTD*.all.parquet** file.
* **QTD\*.permuted.parquet** - Permutation p-values for all tested molecular traits.
* **QTD\*.credible_sets.parquet** - Fine mapped credible sets detected by SuSiE (purity filtered).
* **QTD\*.lbf_variable.parquet** - Fine mapped log Bayes factors from SuSiE for all conditinally distinct signals.


# Column names of the nominal cis-QTL summary statistics files (QTD*.all.parquet and QTD*.cc.parquet)

* **molecular_trait_id** - ID of the molecular trait used for QTL mapping. Depending on the quantification method used, this can be either a gene id, exon id, transcript id or a txrevise promoter, splicing or 3'end event id. Examples: ENST00000356937, ENSG00000008128.
* **chromosome** - GRCh38 chromosome name of the variant (e.g. 1,2,3 ...,X).
* **position** - GRCh38 position of the variant.
* **ref** - GRCh38 reference allele.
* **alt** - GRCh38 alternative allele (also the effect allele).
* **variant** - The variant ID (chromosome_position_ref_alt) e.g. chr19_226776_C_T. Based on GRCh38 coordinates and reference genome. The chromosome, position, ref and alt values should exactly match same fields in the summary statistics file, with 'chr' prefix added to the chromosome number. 
* **ma_samples** - Number of samples carrying at least one copy of the minor allele.
* **maf** - Minor allele frequency within a sample group (e.g. cell type or tissues within a study).
* **pvalue** - Nominal p-value of association between the variant and the molecular trait.
* **beta** - Regression coefficient from the linear model.
* **se** - Standard error of the beta.
* **type** - Type of the genetic variant (SNP, INDEL or OTHER).
* **ac** - Count of the alternative allele.
* **an** - Total number of alleles. For autosomal variants, this is usually two times the sample size. Conversly, for autosomal variants, sample size is equal to an/2.
* **r2** - Optional imputation quality score from the imputation software, can be replaced with NA if not available.
* **molecular_trait_object_id** - For phenotypes with multiple correlated alternatives (multiple alternative transcripts or exons within a gene, multple alternative promoters in txrevise, multiple alternative intons in Leafcutter), this defines the level at which the phenotypes were aggregated. Permutation p-values are calculated accross this set of alternatives.  
* **gene_id** - Ensembl gene ID of the molecular trait. 
* **median_tpm** - Median transcripts per million (TPM) expression value of the gene. Can be replaced with NA if not availble (e.g. in microarray studies).
* **rsid** - The dbSNP v151 rsid of the variant. If the same variant has multiple rsids then these should be split over multiple rows so that all of the other values are duplicated.

# Column names of the permutation p-value files (QTD*.permuted.parquet)

* **molecular_trait_object_id** - For phenotypes with multiple correlated alternatives (multiple alternative transcripts or exons within a gene, multple alternative promoters in txrevise, multiple alternative intons in Leafcutter), this defines the level at which the phenotypes were aggregated. Permutation p-values are calculated accross this set of alternatives.
* **molecular_trait_id** - ID of the molecular trait used for QTL mapping. Depending on the quantification method used, this can be either a gene id, exon id, transcript id or a txrevise promoter, splicing or 3'end event id. Examples: ENST00000356937, ENSG00000008128. 
* **n_traits** - The number of molecular traits over which permutation p-values were calculated (e.g. the number of transcripts per gene). Note that the permutations are performed accross all molecular traits within the same molecular trait object (e.g. all transcripts of a gene) and the results are reported for the most significant variant and molecular trait pair. 
* **n_variants** - number of genetic variants tested within the cis region of the molecular trait.
* **variant** - The variant ID (chromosome_position_ref_alt) e.g. chr19_226776_C_T. Based on GRCh38 coordinates and reference genome. The chromosome, position, ref and alt values should exactly match same fields in the summary statistics file, with 'chr' prefix added to the chromosome number.
* **chromosome** - GRCh38 chromosome name of the variant (e.g. 1,2,3 ...,X).
* **position** - GRCh38 position of the variant.
* **pvalue** - Nominal p-value of association between the variant and the molecular trait.
* **beta** - Regression coefficient from the linear model.
* **p_perm** - Empirical p-value calculated from 1000 permutations.
* **p_beta** - Estimated empirical p-value based on the beta distribution. This is the column that you want to use for filtering the results. See the FastQTL [paper](http://dx.doi.org/10.1093/bioinformatics/btv722) for more details. 

# Column names of the fine mapped credible set files from SuSiE (QTD*.credible_sets.parquet)

* **molecular_trait_id** - ID of the molecular trait used for QTL mapping. Depending on the quantification method used, this can be either a gene id, exon id, transcript id or a txrevise promoter, splicing or 3'end event id. Examples: ENST00000356937, ENSG00000008128.
* **chromosome** - GRCh38 chromosome name of the variant (e.g. 1,2,3 ...,X).
* **position** - GRCh38 position of the variant.
* **ref** - GRCh38 reference allele.
* **alt** - GRCh38 alternative allele (also the effect allele).
* **variant** - The variant ID (chromosome_position_ref_alt) e.g. chr19_226776_C_T. Based on GRCh38 coordinates and reference genome. The chromosome, position, ref and alt values should exactly match same fields in the summary statistics file, with 'chr' prefix added to the chromosome number.  
* **ma_samples** - Number of samples carrying at least one copy of the minor allele.
* **maf** - Minor allele frequency within a sample group (e.g. cell type or tissues within a study).
* **pvalue** - Marginal nominal p-value of association between the variant and the molecular trait.
* **beta** - Marginal regression coefficient from the linear model.
* **se** - Standard error of the beta.
* **type** - Type of the genetic variant (SNP, INDEL or OTHER).
* **ac** - Count of the alternative allele.
* **r2** - Optional imputation quality score from the imputation software, can be replaced with NA if not available.
* **molecular_trait_object_id** - For phenotypes with multiple correlated alternatives (multiple alternative transcripts or exons within a gene, multple alternative promoters in txrevise, multiple alternative intons in Leafcutter), this defines the level at which the phenotypes were aggregated. Permutation p-values are calculated accross this set of alternatives.
* **gene_id** - Ensembl gene ID of the molecular trait. 
* **median_tpm** - Median transcripts per million (TPM) expression value of the gene. Can be replaced with NA if not availble (e.g. in microarray studies).
* **rsid** - The dbSNP v151 rsid of the variant. If the same variant has multiple rsids then these should be split over multiple rows so that all of the other values are duplicated.
* **cs_id** - Unique ID for each credible set within a dataset.
* **cs_size** - Size of the credible set.
* **pip** - Posterior inclusion probability the variant.
* **z** - Univariate z-score for the variant (calculated by SuSiE).
* **cs_min_r2** - minimal LD (r2) between any two variants in the credible set
* **region** - start and end coordinates of the fine mapped region
