library("data.table")
library("dplyr")
library("devtools")
devtools::load_all("../eQTLUtils/")

#' Import transcript metadata from biomart web export
#'
#' @param biomart_path Path to the biomart download text file
#' @param col_types Column types of the biomart download file (for readr)
#'
#' @return tibble containing transcript metadata
#' @export
importBiomartMetadata <- function(biomart_path, col_types = "ccccccccciiciiicccccccccidccccii"){
  transcript_meta = readr::read_tsv(biomart_path, col_types = col_types)
  col_df = dplyr::data_frame(column_name = c('Gene stable ID', 'Transcript stable ID', 'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)', 'Strand', 'Transcript start (bp)', 'Transcript end (bp)', 'Transcription start site (TSS)', 'Transcript length (including UTRs and CDS)', 'Transcript support level (TSL)', 'APPRIS annotation', 'GENCODE basic annotation', 'Gene name', 'Transcript name', 'Transcript count', 'Transcript type', 'Gene type', 'Gene % GC content', 'Version (gene)', 'Version (transcript)'),
                             column_id = c('gene_id', 'transcript_id', 'chromosome', 'gene_start', 'gene_end', 'strand', 'transcript_start', 'transcript_end', 'tss', 'transcript_length', 'transcript_tsl', 'transcript_appris', 'is_gencode_basic', 'gene_name', 'transcript_name', 'transcript_count', 'transcript_type', 'gene_type', 'gene_gc_content', 'gene_version', 'transcript_version'))
  transcript_meta = transcript_meta[,col_df$column_name] %>% dplyr::distinct()
  colnames(transcript_meta) = col_df$column_id
  transcript_meta = dplyr::mutate(transcript_meta, transcript_tsl = ifelse(transcript_tsl == "tslNA", NA, transcript_tsl))
    
  return(transcript_meta)
}

extractGeneMetadataFromBiomartFile <- function(biomart_df){
  required_gene_meta_columns = c("phenotype_id","quant_id","group_id","gene_id",
                                 "chromosome","gene_start","gene_end","strand",
                                 "gene_name","gene_type","gene_gc_content",
                                 "gene_version","phenotype_pos")
  
  #Extract gene metadata
  gene_data = dplyr::select(biomart_df, gene_id, chromosome, gene_start, gene_end, strand,
                            gene_name, gene_type, gene_gc_content, gene_version) %>%
    dplyr::distinct() %>%
    dplyr::mutate(phenotype_id = gene_id, group_id = gene_id, quant_id = gene_id) %>%
    dplyr::mutate(phenotype_pos = ifelse(strand == 1, gene_start, gene_end)) %>%
    dplyr::select(all_of(required_gene_meta_columns), dplyr::everything())
  return(gene_data)
}

#### Gene counts (featureCounts) ####
#Specify required phenotype metadata columns
required_phenotype_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                               "gene_end","strand","gene_name","gene_type","gene_version","phenotype_pos")
required_gene_meta_columns = c(required_phenotype_meta_columns, "phenotype_gc_content", "phenotype_length")

#Import gene metadata
transcript_meta = importBiomartMetadata("data/Homo_sapiens.GRCh38.105_biomart_download.txt.gz")
gene_metadata = extractGeneMetadataFromBiomartFile(transcript_meta) %>%
  dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

#Import featureCounts gene lengths
lengths = read.table("data/aipt_A.sorted_gene.featureCounts.txt", 
                     header = T, stringsAsFactors = F) %>% 
  dplyr::as_tibble() %>%
  dplyr::transmute(phenotype_id = Geneid, phenotype_length = Length) %>%
  dplyr::filter(!(phenotype_id %like% "PAR_Y")) %>%
  reformatPhenotypeId()

#Make gene metadata
final_gene_metadata = dplyr::left_join(lengths, gene_metadata, by = "phenotype_id") %>%
  dplyr::rename(phenotype_gc_content = gene_gc_content) %>%
  dplyr::select(required_gene_meta_columns, dplyr::everything()) 

#Save expression matrix
gz2 = gzfile("data/phenotype_metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz", "w")
write.table(final_gene_metadata, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Exon counts ####
exon_nuc_content = readr::read_tsv("data/DEXSeq/gencode.v39.annotation.nochr.patched_contigs.DEXSeq.nuc_content.txt",
                                   col_types = "ccciiccccddiiiiiii")
nuc_content = exon_nuc_content[,c(1,4,5,11)]
colnames(nuc_content) = c("seqnames", "start", "end", "phenotype_gc_content")
nuc_content = dplyr::distinct(nuc_content)

#Subset of gene metadata
gene_metadata_subset = dplyr::select(gene_metadata, gene_id, chromosome, gene_start, gene_end, 
                                     strand, gene_name, gene_type, gene_version)
  
#Make exon metadata
exon_metadata = rtracklayer::import.gff("data/DEXSeq/gencode.v39.annotation.nochr.patched_contigs.DEXSeq.gff") %>%
  as.data.frame() %>% dplyr::as_tibble() %>%
  dplyr::filter(type != "aggregate_gene") %>%
  dplyr::mutate(phenotype_id = paste(gene_id, seqnames, start, end, sep = "_")) %>%
  dplyr::mutate(phenotype_pos = floor((start + end)/2)) %>%
  dplyr::transmute(phenotype_id, gene_id, seqnames, start, end, phenotype_length = width, phenotype_pos) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  dplyr::left_join(nuc_content, by = c("seqnames", "start", "end")) %>%
  dplyr::filter(!grepl("+", gene_id, fixed = TRUE)) %>%
  removeGeneVersion() %>%
  dplyr::mutate(quant_id = gene_id, group_id = gene_id) %>%
  dplyr::left_join(gene_metadata_subset, by = "gene_id") %>%
  dplyr::select(-seqnames, -start, -end) %>%
  dplyr::select(required_gene_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("data/phenotype_metadata/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz", "w")
write.table(exon_metadata, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)


#### GENCODE transcripts ####
tx_gene_map = dplyr::transmute(transcript_meta, phenotype_id = transcript_id, gene_id)
tx_estimates = readr::read_tsv("data/gencode.v39.transcripts.TPM.merged.txt") %>% reformatPhenotypeId()
gencode_transcripts_meta = dplyr::select(tx_estimates, phenotype_id) %>% 
  dplyr::left_join(tx_gene_map) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("data/phenotype_metadata/transcript_usage_Ensembl_105_phenotype_metadata.tsv.gz", "w")
write.table(gencode_transcripts_meta, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)


### txrevise events ####
list = as.list(c("data/txrevise.grp_1.upstream.gff3", 
                 "data/txrevise.grp_2.upstream.gff3",
                 "data/txrevise.grp_1.contained.gff3",
                 "data/txrevise.grp_2.contained.gff3",
                 "data/txrevise.grp_1.downstream.gff3",
                 "data/txrevise.grp_2.downstream.gff3"))
event_quants = purrr::map_df(list, ~rtracklayer::import.gff3(.) %>% 
                               dplyr::as_tibble() %>% 
                               dplyr::filter(type == "mRNA"))
event_ids = dplyr::transmute(event_quants, phenotype_id = ID)

#Extract relevant infromation from the txrevise events
txrevise_meta = dplyr::select(event_ids, phenotype_id) %>% 
  tidyr::separate(phenotype_id, c("gene_id", "grp", "pos", "transcript"), sep = "\\.", remove = F) %>%
  dplyr::mutate(quant_id = paste(gene_id, grp, pos, sep = "."), group_id = paste(gene_id, pos, sep = ".")) %>% 
  dplyr::select(phenotype_id, quant_id, group_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content, -group_id, -quant_id), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("data/phenotype_metadata/txrevise_Ensembl_105_phenotype_metadata.tsv.gz", "w")
write.table(txrevise_meta, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Human HT-12 V4 metadata ####
human_ht12 = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz") %>%
  dplyr::select(phenotype_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything()) %>%
  dplyr::filter(!is.na(gene_name))

#Save metadata file
gz2 = gzfile("metadata/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(human_ht12, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Human HT-12 V4 metadata ####
affy = readr::read_tsv("metadata/gene_metadata/Affy_Human_Gene_1_0_ST_gene_metadata.txt.gz", col_types = "ccccciiiccddi") %>%
  dplyr::select(phenotype_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything()) %>%
  dplyr::filter(!is.na(gene_name))

#Save metadata file
gz2 = gzfile("metadata/phenotype_metadata/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(affy, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)



