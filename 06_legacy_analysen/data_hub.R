####################
# central data hub #
####################
# used to merge and match biobank samples throughout the different runs

# whole GMALL biobank registry
whole_biobank <- readxl::read_excel("Data/GMALL_biobank_probenuebersicht_b_precursor.xlsx", sheet = 2)

# samples from the n=100 pilot run
samples_pilotrun <- as.character(meta.all$ID)

# samples from the n=79 main run
samples_mainrun <- c("76", "1773", "552", "1769", "2293", "3273", "3054", "1647", "2118", "3437", "517", "1912", "15-0456-01", "1043", "3379", "1025", "3209", "600", "1935", "2149", "1743", "1619", "3684", "3478", "2879", "2178", "1491", "1297", "1805", "2134", "1007", "1660", "1656", "562", "3366", "690", "3831", "1585", "2361", "1756", "272", "962", "1782", "3987", "2190", "2021", "3022", "1554", "2314", "1559", "4842", "892", "1092", "3147", "4539", "2511", "144", "2893", "3840", "261", "1524", "1035", "837", "473", "155", "2962", "625", "372", "1866", "1738", "2530", "862", "502", "442", "4579", "3978", "3581", "476", "4469", "859", "4847", "531", "4535", "51", "3632", "1122", "3264", "4945", "2346", "4871", "3476", "2504", "3585", "843", "1024", "3949", "3506", "3652", "4759", "2406", "4416", "3744", "2477", "538", "2492", "4382", "3497", "3226", "4776", "1480", "1801", "1801", "3607", "3020", "4658", "4493", "4553", "320", "4020", "1183", "3630", "1183", "740", "2149", "874", "760", "1079", "2178", "699", "559", "3524", "1277", "1277", "4071", "4826", "4167", "3762", "3756", "3722", "3720", "3400", "3312", "2873", "2852", "2703", "2700", "2688", "2687", "2630", "2599", "2570", "2512", "1944", "1385")


# label biobank samples that were already measured in any of the two runs

whole_biobank %>% 
  mutate(measured_in_main = factor(BM_Nr %in% samples_mainrun)) %>% 
  relocate(c(measured_in_main), .after = ID) -> whole_biobank

table(whole_biobank$measured_in_main)

whole_biobank %>% 
  mutate(across(.cols = c(N_Prob_KM, N_Prob_PB), ~ ifelse(is.na(.), 1, .))) -> whole_biobank

whole_biobank %>% 
  dplyr::select(ID, measured_in_main) %>% 
  clipr::write_clip()


#### 2024/03/26
# match with samples measured in HD 

hd_ids <- readxl::read_excel("~/Forschung/AG Oellerich/ALL Proteomics/Data/Raffel_Rhapsody_GMALL_identifiers.xlsx")
hd_ids %>% 
  filter(cohort == "F") -> hd_ids

# try matching using the lab ID
# sample_metadata is the metadata frame from the 2nd cohort
# hd uses the BM_Nr as the sample_name_corrected

ggvenn::ggvenn(data = list("FFM_MS_cohort" = sample_metadata$lab_id, "HD_Rhapsody_cohort" = hd_ids$LabID))

full_join(hd_ids, dplyr::select(sample_metadata, ms_name, sample_id, sample_name, ID, lab_id), by=c("LabID" = "lab_id"), multiple = "first")

ggvenn::ggvenn(data = list("FFM_MS_cohortI" = sample_metadata$sample_name, "HD_Rhapsody_cohort" = hd_ids$sample_name_corrected))

# there are 57 samples that don't map

full_join(hd_ids, dplyr::select(sample_metadata, ms_name, sample_id, sample_name, ID, lab_id), by=c("LabID" = "lab_id"), multiple = "first") %>% 
  filter(is.na(rhapsody_id) |is.na(sample_name)) %>% 
  left_join(dplyr::select(whole_biobank, ID, BM_Nr, measured_in_main) %>% dplyr::rename("biobank_ID" = "ID"), by=c("sample_name_corrected" = "BM_Nr"), multiple = "first") %>% 
  left_join(dplyr::select(sample_metadata, ms_name, sample_id, ID), by=c("biobank_ID" = "ID"), multiple = "first")


# we need to identify 19 samples that we measured but are not in the list

filter(sample_metadata, !lab_id %in% hd_ids$LabID) -> ffm_samples_not_in_hd

# hd F-assigned samples not in FFM
filter(hd_ids, cohort == "F" & !LabID %in% sample_metadata$lab_id) -> hd_samples_not_in_ffm

hd_samples_not_in_ffm$sample_name_corrected %in% sample_metadata

# comparison of masterfile to biobank overview

masterfile <- readxl::read_excel("Data/Rhapsody_GMALL_identifiers 25MAR2024_sent to FFM_AC.xlsx", sheet = "FFM_Samples_not_in_HD_cohort")

masterfile %>% 
  left_join(sample_metadata, by=c("Lab_ID" = "lab_id"), na_matches = "never") %>% 
  mutate(ms_cohortII = !is.na(ms_name)) %>% 
  dplyr::select(-ID) %>% 
  mutate(Raffel = Lab_ID %in% masterfile_raffel$LabID) %>% 
  clipr::write_clip()

#check that no labds from cohortII are missing
filter(sample_metadata, !lab_id %in% masterfile$Lab_ID)

# raffel-ids

masterfile_raffel <- readxl::read_excel("Data/Rhapsody_GMALL_identifiers 25MAR2024_sent to FFM_AC.xlsx", sheet = 1)

sum(masterfile_raffel$LabID %in% masterfile$Lab_ID)

mutate(masterfile_raffel, in_ffm = !LabID %in% masterfile$Lab_ID & !is.na(LabID)) %>% 
  clipr::write_clip()

bm_nr_to_labid <- readxl::read_excel("Data/Metadata/Kiel_RNASeq_molecular_classes_all_mainrun.xlsx")

# masterfile anjali

masterfile_ac_cohortI <- readxl::read_excel("Data/Masterfile Samples B-ALL Kohorten_SW_AC.xlsx", sheet = 1, range = c("A1:M14"))
masterfile_ac_cohortII <- readxl::read_excel("Data/Masterfile Samples B-ALL Kohorten_SW_AC.xlsx", sheet = 2, range = c("A1:V87")) 
masterfile_ac_cohortIII <- readxl::read_excel("Data/Masterfile Samples B-ALL Kohorten_SW_AC.xlsx", sheet = 3, range = c("A1:T68")) 

raffel_ids <- c("445nkg", "383gda", "741ojf", "097sle", "368ogl", "223vem", "184boa", "147pqr", "369spm", "604kos", "059kot", "261esh", "791dxf", "396elr", "432eff", "372ebl", "616doz", "260jzp", "790cmw", "802kue", "876xud", "630pyv", "937fms", "358huu", "160thi", "532jxa", "075uuh", "581eib", "025dtk", "357cgw", "630hlq", "876dcp", "616bzy", "628ogm", "962lqj", "050ctw", "109spm", "111uqm", "677est", "780nho", "925xel", "582bpc", "306kqj", "505ivi", "233dpb", "974lej", "840cjr", "394ncf", "949bcd", "074jih", "394cdf", "147wyz", "162wtd", "655oyt", "715swa", "727qoz", "642wae", "594ajd", "369uxn", "346dax", "925xex", "432twn", "901qti", "208qwq", "628jhs", "310hzs", "628jtd", "074ejo", "815ldh", "738hah", "308kzp", "482iiw", "467cmv", "579pjc", "666kvh", "627pny", "108oxr", "013jlm", "741boo", "924mvr", "902fhi", "678dqd", "051fxp", "085bxt", "665hbu", "531ukm", "418zej", "161edq", "679qbn", "631maa", "471smg", "689iqx", "505ebm", "741hec", "776blu", "838poy", "219fai", "754xla", "135jdx", "147lkj", "247mdt", "752clj", "134ymo", "049dzk", "209qwq", "418bdq", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

tibble(lab_id = raffel_ids) %>% 
  left_join(masterfile_ac_cohortII, by=c("lab_id" = "LabID"), na_matches = "never") %>% 
  clipr::write_clip()

### creating a rna overview for the Kiel sample sheet from all available samples
### goal is one unique rna sample per labID 

masterfile_ac_cohortI %>% 
  dplyr::select(`BM Nr`, LabID, RNA, RIN...11, `RNA conc (ng/ul)`) %>% 
  mutate(cohort = "cohort I") -> masterfile_ac_cohortI_RNA

masterfile_ac_cohortII %>% 
  dplyr::select(`BM Nr`, LabID, RNA, RIN, `RNA conc (ng/ul)`) %>% 
  mutate(cohort = "cohort II",
         `BM Nr` = as.numeric(`BM Nr`)) -> masterfile_ac_cohortII_RNA

masterfile_ac_cohortIII %>% 
  dplyr::select(`BM Nr`, LabID, RNA, RIN, `RNA conc (ng/ul)`) %>% 
  mutate(cohort = "cohort III") -> masterfile_ac_cohortIII_RNA

bind_rows(masterfile_ac_cohortI_RNA, masterfile_ac_cohortII_RNA, masterfile_ac_cohortIII_RNA) %>%
  mutate(LabID = ifelse(LabID == "na" & cohort == "cohort I", paste0("cohortI_bmnr_", `BM Nr`), LabID),
         LabID = ifelse(`BM Nr` %in% c(1866L, 1297L, 2295L, 3980L), paste0("bmnr_", `BM Nr`), LabID)) -> masterfile_ac_cohort_RNA

ids_gesamtübersicht <- c("025dtk", "025dtk", "049dzk", "050ctw", "050ctw", "050ctw", "051fxp", "059kot", "074ejo", "074jih", "075uuh", "075uuh", "085bxt", "097sle", "097sle", "108oxr", "109spm", "111uqm", "134ymo", "135jdx", "147lkj", "147pqr", "147wyz", "160thi", "161edq", "162wtd", "184boa", "184boa", "208qwq", "209qwq", "219fai", "223vem", "233dpb", "233dpb", "247mdt", "260jzp", "260jzp", "261esh", "306kqj", "308kzp", "310hzs", "310hzs", "346dax", "357cgw", "358huu", "368ogl", "368ogl", "369spm", "369uxn", "372ebl", "383gda", "383gda", "383gda", "383gda", "394cdf", "394ncf", "396elr", "396elr", "418bdq", "418zej", "418zej", "432eff", "432eff", "432eff", "432eff", "432twn", "432twn", "432twn", "445nkg", "467cmv", "471smg", "482iiw", "505ebm", "505ivi", "531ukm", "532jxa", "532jxa", "579pjc", "581eib", "582bpc", "582bpc", "582bpc", "594ajd", "604kos", "604kos", "604kos", "616bzy", "616doz", "616doz", "628jhs", "628jtd", "628ogm", "628ogm", "630hlq", "630hlq", "630pyv", "630pyv", "631maa", "655oyt", "665hbu", "666kvh", "678dqd", "679qbn", "689iqx", "715swa", "727qoz", "738hah", "738hah", "741boo", "741hec", "741ojf", "741ojf", "752clj", "754xla", "776blu", "790cmw", "790cmw", "791dxf", "802kue", "815ldh", "838poy", "840cjr", "840cjr", "876dcp", "876xud", "901qti", "901qti", "902fhi", "924mvr", "925xel", "925xex", "937fms", "949bcd", "962lqj", "962lqj", "974lej")


masterfile_ac_cohort_RNA %>% 
  mutate(duplette = LabID %in% duplicated(masterfile_ac_cohort_RNA$LabID)) %>% 
  mutate(duplette = ifelse(duplette, "Duplette", NA_character_)) %>% 
  filter(LabID %in% unique(ids_gesamtübersicht) | stringr::str_detect(LabID, "cohort") |  `BM Nr` %in% c(1866L, 1297L, 2295L, 3980L)) %>% 
  group_by(LabID) %>% 
  mutate(index = row_number()) %>% 
  arrange(LabID) %>% 
  pivot_wider(id_cols = LabID, names_from = index, values_from = c(`RNA conc (ng/ul)`, RIN)) %>% 
  rename_with(.cols = everything(), .fn = ~ stringr::str_replace_all(., " conc \\(ng/ul\\)_", "_conc_")) %>% 
  dplyr::relocate(RIN_1,.after = RNA_conc_1) %>% 
  dplyr::relocate(RIN_2, .after = RNA_conc_2) %>% 
  dplyr::relocate(RIN_3, .after = RNA_conc_3) %>% 
  
  

