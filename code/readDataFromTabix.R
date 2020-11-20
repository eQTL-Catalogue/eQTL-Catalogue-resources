library("dplyr")

import_eQTLCatalogue <- function(ftp_path, region, column_names, molecular_trait_id, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with Rsamtools
  mtid = molecular_trait_id
  summary_stats = scanTabixDataFrame(ftp_path, region, col_names = column_names)[[1]]

  #Remove rsid duplicates and multi-allelic variant
  if(!is.null(summary_stats)){
    summary_stats = dplyr::filter(summary_stats, molecular_trait_id == mtid) %>%
      dplyr::select(-rsid) %>% 
      dplyr::distinct() %>% #rsid duplicates
      dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
      dplyr::group_by(id) %>% 
      dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
      dplyr::filter(row_count == 1) #Multialllics 
  } else{
    summary_stats = NULL
  }
  
  return(summary_stats)
}

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

importTraitVariantPair <- function(molecular_trait_id, variant, ftp_path){
  
  #Read column names
  column_names = colnames(read.table(ftp_path, nrows = 1, header = T))
  
  #Extract chromosme and position from variant id
  variant_split = stringr::str_split(variant, pattern = "_") %>% unlist()
  chromosome = stringr::str_replace(variant_split[1], "chr", "")
  position = as.integer(variant_split[2])
  
  #Make region GRanges object
  region_granges = GenomicRanges::GRanges(
    seqnames = chromosome, 
    ranges = IRanges::IRanges(start = position, end = position), 
    strand = "*")
  
  result = import_eQTLCatalogue(ftp_path, region_granges, column_names, molecular_trait_id)
  return(result)
}

#Find all files
ftp_list = list.files("/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/", pattern = "*exon.nominal.sorted.tsv.gz$")
file_list = setNames(as.list(paste0("/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/", ftp_list)), ftp_list)

#Extract effects
results_df = purrr::map_df(file_list, ~importTraitVariantPair("ENSG00000113161.16_5_75355215_75355364", "chr5_75355259_A_G", .), .id = "file_name")
write.table(results_df, "HMGCR_exon.tsv", sep = "\t", row.names = F, quote = F)

#Find all files
ftp_list = list.files("/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/", pattern = "*ge.nominal.sorted.tsv.gz$")
file_list = setNames(as.list(paste0("/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/", ftp_list)), ftp_list)

#Extract effects
results_df = purrr::map_df(file_list, ~importTraitVariantPair("ENSG00000113161", "chr5_75355259_A_G", .), .id = "file_name")
write.table(results_df, "HMGCR_gene.tsv", sep = "\t", row.names = F, quote = F)


