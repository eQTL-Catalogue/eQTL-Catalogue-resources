library("readr")
library("dplyr")
library("rtracklayer")
library("UpSetR")

read_colocs_to_df <- function(file_dir, ...) {
  file_list <- list.files(file_dir, recursive = TRUE ,full.names = TRUE)
  dataset = data.frame()
  
  for (file in file_list){
    print(file)    
    
    # if the merged dataset does exist, append to it
    temp_dataset <-read_tsv(file, trim_ws = TRUE)
    temp_dataset <- temp_dataset %>% dplyr::filter(...)
    dataset<-rbind(dataset, temp_dataset)
  }
  return(dataset)
}

#Import molecular trait metadata
tx_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::transmute(molecular_trait_id = phenotype_id, gene_id)
exon_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/exon_counts_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::transmute(molecular_trait_id = phenotype_id, gene_id)
txrev_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/txrevise_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::transmute(molecular_trait_id = phenotype_id, gene_id)

#Import all coloc
ge_colocs = read_colocs_to_df("coloc/ge/") %>%
  dplyr::filter(PP.H4.abf > 0.5) %>%
  dplyr::mutate(gene_id = molecular_trait_id) %>%
  dplyr::mutate(quant = "ge")
tx_colocs = read_colocs_to_df("coloc/tx/") %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::left_join(tx_meta, by = "molecular_trait_id") %>%
  dplyr::mutate(quant = "tx")
exon_colocs = read_colocs_to_df("coloc/exon/") %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::left_join(exon_meta, by = "molecular_trait_id") %>%
  dplyr::mutate(quant = "exon")
txrev_colocs = read_colocs_to_df("coloc/txrev/") %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::left_join(txrev_meta, by = "molecular_trait_id") %>%
  dplyr::mutate(quant = "txrev")
joint_colocs = dplyr::bind_rows(ge_colocs, tx_colocs, exon_colocs, txrev_colocs)

#Count HMGCR colocs
dplyr::filter(tx_colocs, gene_id == "ENSG00000113161")
dplyr::filter(exon_colocs, gene_id == "ENSG00000113161")
dplyr::filter(txrev_colocs, gene_id == "ENSG00000113161")

#Import inependet LD blocks
regions = readr::read_tsv("../novel_coloc_plots/ind_regions_Hg38.bed")
regions_ranges = GenomicRanges::GRanges(seqnames = as.character(regions$chr), ranges = IRanges::IRanges(start = regions$start, end = regions$end))
elementMetadata(regions_ranges) = regions

#Make variant GRanges object
variant_regions = joint_colocs %>%
  dplyr::select(variant) %>% 
  dplyr::distinct() %>% 
  tidyr::separate(variant, c("chr", "pos", "ref", "alt"), sep = "_", remove = F) %>% 
  dplyr::mutate(chr = stringr::str_replace(chr, "chr", "")) %>% 
  dplyr::mutate(pos = as.integer(pos)) %>%
  dplyr::arrange(chr, pos)
variant_ranges = GenomicRanges::GRanges(seqnames = as.character(variant_regions$chr), ranges = IRanges::IRanges(start = variant_regions$pos, end = variant_regions$pos))
elementMetadata(variant_ranges) = variant_regions

#Find overlaps
olaps = findOverlaps(variant_ranges, regions_ranges)
olap_variants = variant_regions[queryHits(olaps),]
olap_variants$ld_block_id = elementMetadata(regions_ranges[subjectHits(olaps),])$ID_hg38
olap_variants = dplyr::select(olap_variants, variant, ld_block_id)

#Count the number of colocs for each quantification method
counts = dplyr::left_join(joint_colocs, olap_variants, by = "variant") %>%
  dplyr::transmute(ld_block_id, quant, count = 1) %>% 
  dplyr::distinct() %>% 
  tidyr::pivot_wider(names_from = quant, values_from = count) %>%
  dplyr::mutate(ge = ifelse(is.na(ge), 0, ge)) %>%
  dplyr::mutate(tx = ifelse(is.na(tx), 0, tx)) %>%
  dplyr::mutate(exon = ifelse(is.na(exon), 0, exon)) %>%
  dplyr::mutate(txrev = ifelse(is.na(txrev), 0, txrev)) %>%
  dplyr::rename(gene = ge, txrevise = txrev)



#Make UpSetR plot
pdf(file="LDLC_upet_plot.pdf", width = 3.5, height = 2.3, onefile=FALSE)
upset(as.data.frame(counts), sets = c("gene","tx","exon", "txrevise"), order.by = "freq", mainbar.y.label	= "Number of colocalisations")
dev.off()
