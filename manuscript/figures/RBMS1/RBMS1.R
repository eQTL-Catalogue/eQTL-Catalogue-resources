library("dplyr")
library("gwasvcf")
library("ggplot2")

#Visualise effect sizes
#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param An instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with Rsamtools
  summary_stats = scanTabixDataFrame(ftp_path, region, col_names = column_names)[[1]] %>%
    dplyr::filter(gene_id == selected_gene_id)
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::select(summary_stats, -rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
  
  return(summary_stats)
}



#Import summary stats
eqtl_data = readr::read_tsv("BLUEPRINT_PE.T-cell_ge.nominal.sorted.tsv.gz") %>%
  dplyr::filter(molecular_trait_id == "ENSG00000153250") %>%
  dplyr::mutate(LP = -log(pvalue,10)) %>%
  dplyr::select(position, LP) %>%
  dplyr::mutate(type = "RBMS1 eQTL")
gwas_data = gwasvcf::query_gwas("LC_GWAS_subset.GRCh38.sorted.vcf.gz", chrompos = "2:160100000-160700000") %>% 
  gwasvcf::vcf_to_granges() %>% 
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(position = start, LP, type = "Lymphocyte count")
gwas_data_selected = dplyr::semi_join(gwas_data, eqtl_data, by = "position")

joint_data = dplyr::bind_rows(eqtl_data, gwas_data)

#Import eQTL credible sets
#cs = readr::read_tsv("BLUEPRINT_PE.T-cell_ge.purity_filtered.txt.gz") %>%
#  dplyr::filter(phenotype_id == "ENSG00000153250")
#write.table(cs, "RBMS1_credible_set.tsv", sep = "\t", quote = F, row.names = F)
cs = readr::read_tsv("RBMS1_credible_set.tsv")

#Flag the credible set in GWAS results
gwas_flagged = dplyr::mutate(joint_data, in_cs = ifelse(position %in% cs$pos, TRUE, FALSE))

#Make manhattan plots
manhattan = ggplot(gwas_flagged, aes(y = LP, x = position, color = in_cs)) + 
  geom_point() + 
  facet_grid(type~.) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Chromosome 2 position") + 
  ylab("-log10(p-value)") + 
  scale_colour_manual(name = "group",
                      values=c("black","red")) +
  theme(legend.position = "none")
ggsave("RBMS1_manhattan.pdf", plot = manhattan, width = 4, height = 3)



# Make colocalisation figure
#Import all coloc results
#files_list = list.files("LC-ebi-a-GCST004627_rnaseq/")
#paths = paste0("LC-ebi-a-GCST004627_rnaseq/", files_list)
#path_list = setNames(as.list(paths), files_list)
#coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
#rbms1_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T") %>%
#  dplyr::mutate(qtl_subset = stringr::str_replace(qtl_subset, pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT")) %>%
#  dplyr::mutate(qtl_subset = stringr::str_replace(qtl_subset, pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT"))


#Import all coloc results
#files_list = list.files("LC-ebi-a-GCST004627_GTEx_V8/")
#paths = paste0("LC-ebi-a-GCST004627_GTEx_V8/", files_list)
#path_list = setNames(as.list(paths), files_list)
#coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
#gtex_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T")

#rbms1_colocs = dplyr::bind_rows(rbms1_colocs, gtex_colocs)
#write.table(rbms1_colocs, "rbms1_coloc_results.tsv", sep = "\t", row.names = F, quote = F)
rbms1_colocs = read.table("rbms1_coloc_results.tsv", header = TRUE, stringsAsFactors = F)

#Visualise
coloc_plot = ggplot(rbms1_colocs, aes(x = PP.H3.abf, y = PP.H4.abf, label = qtl_subset)) + 
  geom_point() + 
  theme_light() +
  theme(panel.grid = element_blank()) +
  xlab("PP3 (distinct causal variants)") +
  ylab("PP4 (shared causal variant)")
ggsave("RBMS1_coloc_plot.pdf", plot = coloc_plot, width = 3, height = 3)

