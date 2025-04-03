# build a final proteogenomic data set
# w/ sample id for which we have genomics
# and append these for which we have proteomics

all_mutations %>% 
  janitor::clean_names() %>% 
  dplyr::select(mll_id, sample_id) %>% 
  dplyr::rename("DNA_mll_id" = "mll_id") %>% 
  mutate(sample_id = stringr::str_to_lower(sample_id)) %>% 
  mutate(omic_DNA = T ) -> av_mutations

rna_counts %>% 
  colnames() %>% 
  stringr::str_remove_all("x") %>% 
  as_tibble() %>% 
  dplyr::rename("sample_id" = "value") %>% 
  mutate(omic_RNA = T) -> av_rna

readxl::read_excel("ALL_samples_sent_for_genomics.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::select(lab_id, bm_nr, id, proteomic_id, cohort) %>% 
  mutate(omic_prot = T) -> all_av_proteomics
  

full_join(av_mutations, av_rna, by=c("sample_id")) -> all_av_genomics

all_av_genomics %>% 
  left_join(all_av_proteomics, by=c("sample_id" = "bm_nr")) -> all_proteogenomics_samplelist

#writexl::write_xlsx(all_proteogenomics_samplelist, "all_proteogenomics_samplelist.xlsx")

