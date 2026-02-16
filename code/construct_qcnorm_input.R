library("dplyr")

dataset_metadata = readr::read_tsv("data_tables/dataset_metadata_upcoming.tsv")

macromap_datasets = dplyr::filter(dataset_metadata, study_id == "QTS000059", quant_method == "ge")
macromap_rnaseq = "/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq/MacroMap"
macromap_sample_meta = "/gpfs/helios/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/MacroMap.tsv"
macromap_vcf = "/gpfs/helios/projects/HipSci/MacroMap/genotypes/MacroMap.vcf.gz"

constructQcnormInputs <- function(dataset_meta, rnaseq_out, sample_meta, vcf){
  inputs = dplyr::select(dataset_meta, study_id, sample_group) %>%
    dplyr::mutate(quant_results_path = file.path(rnaseq_out, sample_group), sample_meta_path = sample_meta, vcf_file = vcf)
  return(inputs)
}

macropmap_inputs = constructQcnormInputs(macromap_datasets, macromap_rnaseq, macromap_sample_meta, macromap_vcf)
write.table(macropmap_inputs, "../SampleArcheology/qcnorm/inputs/MacroMap_inputs.tsv", sep = "\t", row.names = F, quote = F)