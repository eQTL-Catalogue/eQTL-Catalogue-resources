# Metadata standards

## Sample metadata
### Required columns

-   **sample_id** - unique identifier for each sample. Preferably this should the BioSample id (SAMEAXXX) if available. Alternatively run ids from ENA (ERR*) or GEO (GSE*) are also allowed. Finally, if none of the above are available then the original IDs from the study authors are also accepted. 
-   **genotype_id** - unique identifier for each individual in a study. This should exactly match the sample ids in the study VCF file. Note that usually sample_id != genotype_id, because many studies contain multiple samples from the same individual. For RNA-seq data, this should be validated using QTLtools mbv tool, to ensure that no sample swaps have occured.
-   **sex** - Biological sex of the sample. This should be validated based on the gene expression data. Valid values: male, female, NA.
-   **cell_type** - Cell type or tissue of the sample.
-   **condition** - Experimental condition used. For most studies this should be set to ‘naive’.
-  **qtl_group** - Unique combinations of cell types, conditions and other factors. This field is used the to split each dataset into subsets in which each individual occurs only once (i.e. if multiple cell types, stimulation or timepoints are collected from the same set of individuals)
-   **timepoint** - Timepoint of the stimulation, set to 0 if no stimulation was performed.
-   **type** - either RNA-seq or microarray.    
-   **read_length** - RNA-seq read length used (e.g. “75bp”)
-   **stranded** - TRUE if strand-specific RNA-seq protocol was used, FALSE otherwise.    
-   **paired** - TRUE if paired-end sequencing was used, FALSE otherwise.
-  **protocol** - RNA-seq protocol used, valid types are “poly(A)” or “total”. For microarray datasets this should be the microarray platform used (e.g. HumanHT-12_V4 or hugene_10_ST). 
-  **rna_qc_passed** - Sample has passed RNA QC (TRUE/FALSE)
-  **genotype_qc_passed** - Sample has passed genotype QC (TRUE/FALSE) 
-  **study** - unique identifier for the study, usually last name of the first author + year of publication (e.g. Fairfax_2014). The exceptions are studies with established well-known names (GTEx, HipSci, GEUVADIS).
    
### Optional columns
-  **batch** - batch of the RNA-seq or microarray experiment. Used for regressing out batch effects in the microarray datasets.
-  **marker** - Commonly known cell type marker used for sorting cells, such ad CD4 or CD8 for T-cells and CD14 for monocytes.
- **age** - Age of the donor in years.

## Phenotype metadata
### Required columns

 - **phenotype_id** -  id of the molecular trait that has been quantified. This can be either the gene id (RNA-eq eQTLs), probe id (microarray eQTLs), transcript id (full-length transcript usage QTLs), splice junction id (Leafcutter), exon id (exon-level QTLs) or any other molecular trait that has been quantified.
 - **quant_id** - Currently only used for txrevise to normalise promoter/splicing/3'end usage.
-  **group_id** - Used for transcript usage, exon expression and splicing phenotypes. Overlapping phenotypes whose relative expression is quantified belong to the same group (e.g. alternative spliced exons form clusters in Leafcutter).
- **gene_id** - Ensembl gene id. Should be NA if the this information is not available (e.g. novel junction clusters in LeafCutter).
- **gene_start** - Start coordinate of the gene
- **gene_end** - End coordinate of the gene
-   **strand** - Strand of the gene
-   **gene_name** - Gene name extracted from Ensembl biomart.
-   **gene_type** - Gene type (protein coding, lincRNA, etc) extracted from Ensembl biomart.
-  **gene_version** - Ensembl gene version
-  **phenotype_pos** - Genomic position used to determine the centre point of the *cis* window for QTL mapping. 

### Optional columns
- **phenotype_gc_content** - Percentage GC content of the quantified phenotype. Currently used as a covariate by cqn when normalising gene-level and exon-level counts. For genes, this is exracted directly from Ensembl biomart. For exons it is calculated using bedtools nuc command. 
- **phenotype_length** - Length of the phenotype in base pairs. Currently used for gene-level and exon-level counts. Used by cqn and TPM normalisation techniques.

### Definition of the the *cis* window
The center point of the *cis* window is defined by the **phenotype_pos** column in the phenotype metadata file:

 - *gene-level counts* - start of the gene as defined in Ensembl biomart.
 - *exon-level counts* - center point of the exon.
 - *microarray probes* - start of the gene.
 -  *transcript usage* - start of the gene.
 - *txrevise* - start of the gene.
 - leafcutter - center point of the intron cluster.


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTM5MzcyMDE0OSwxMTMzNDkzMzY1LDEzMj
c1NTU5NzNdfQ==
-->