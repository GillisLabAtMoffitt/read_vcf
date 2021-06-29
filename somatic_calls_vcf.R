# Import library
library(tidyverse)
library(data.table)

doParallel::registerDoParallel()
# Load data
# vcf3 <- read_delim(file.choose(), delim = "\t")
data <- readr::read_rds("somatic mutation parsed.rds")

# Data cleaning
somatic_calls <- data %>% 
  mutate_all(funs(str_replace(., ".:.:.:.:.:", NA_character_)))
# somatic_calls1 <- vcf3 %>% 
#   mutate_all(list(str_replace(., ".:.:.:.:.:", NA_character_)))

# somatic_calls1 <- vcf3 %>% 
#   filter_at(vars(starts_with("SL"), str_detect(.,":")))

vcf <- somatic_calls %>% 
  # Unite all mutations characteristics as a row names = make tidy data to allow pivot
  unite(col = "CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT",
        c(1:"FORMAT"), sep = " ") %>% 
  # pivot longer which create 1 row per mutation per patients
  pivot_longer(-CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, 
               names_to = "SLID_tumor", values_to = "DATA") 

# Save each patient ID with or without mutation to add back at the end
unique_patient_in_mutation <- as.data.frame(unique(vcf$patient_id)) %>% 
  `colnames<-` (c("SLID_tumor"))


vcf1 <- vcf %>% 
  # drop the mutations not present (make code faster)
  drop_na("DATA")

saveRDS(vcf1, "mid processed data.rds")

vcf2 <- vcf1 %>% 
  # GT:GQ:GQA:GTS:GTM:DP:AD:DP_Normal:AD_Normal

  # separate back the rownames to have mutation elements
  mutate(mut_element = DATA) %>% 
  separate(mut_element, 
           into = c("GT", "GQ","GQA","GTS","GTM","DP","AD","DP_Normal","AD_Normal"), 
           sep = ":")



saveRDS(vcf2, "int processed data.rds")

vcf3 <- vcf2 %>%
  filter(GQ >= 3)

vcf4 <- vcf3 %>% 
  # separate back the rownames to have mutation elements
  mutate(m_element = CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT) %>% 
  separate(m_element, 
           into = c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), 
           sep = " ") %>% 
  select(-CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, FORMAT)

vcf5 <- vcf4 %>%
  mutate(thousand_genome = str_match(INFO, ";AF_1k=(.*?);")[,2]) %>% # was before "Gene.ensGene=(.*?);GeneDetail.ensGene="
  filter(thousand_genome < 0.01 | is.na(thousand_genome))

saveRDS(vcf5, "int1 processed data.rds")

write_csv(vcf5, "intermediary processed data.csv")


vcf6 <- vcf5 %>%
  mutate(LOCATION = str_match(INFO, ";ANN=(.*?):")[,2]) %>%
  mutate(short = str_match(INFO, ";ANN=(.*?)$")[,2]) %>%
  
  mutate(GENE = str_match(short, ":(.*?)(:|;|,|$)")[,2]) %>% 
  
  mutate(exon = str_match(INFO, "exon(.*?):(.*)")[,2]) %>% 
  
  # VAVIANT_C
  mutate(short = str_match(INFO, "exon(.*?):(.*)")[,3]) %>%
  
  mutate(VARIANT_C = str_match(short, "c.([:alnum:]*)")[,2]) %>%
  # VARIANT_P
  mutate(VARIANT_P = str_match(short, "p.([:alnum:]*)")[,2]) %>%
  
  
  # mutate(gnomAD_noncancer = str_match(INFO, ";non_cancer_AF_popmax=(.*?);")[,2]) %>%
  # # 1. Take the whole string before "esp6500siv2_all"
  # mutate(VARIANT_c = VARIANT_C) %>%
  # # 2. Seaprate by "exon" and coalesce from last to first to get the last "exon" string
  # mutate(VARIANT_c = str_split(string= VARIANT_c, pattern= "exon")) %>%
  # unnest(VARIANT_c) %>%
  # group_by(patient_id, VARIANT_C) %>%
  # mutate(row = row_number()) %>%
  # spread(row, VARIANT_c) %>%
  # ungroup() %>%
  # # separate(VARIANT_c, "exon", into = paste("VARIANT_c", 1:25, sep="_")) %>%
  # mutate(VARIANT_c = coalesce(!!! select(., last_col():"1"))) %>%
  # # 3. Take all between "c". and "p." but will lose data if "p." id not furnished
  # mutate(VARIANT_C = str_match(VARIANT_c, "c.(.*):p.")[,2]) %>%
  # # 4. So do the same with just "c." and coalesce the ones lost (NA) by the "VARIANT_C1" which has them
  # mutate(VARIANT_C1 = str_match(VARIANT_c, "c.(.*)")[,2]) %>%
  # mutate(VARIANT_C = coalesce(VARIANT_C, VARIANT_C1)) %>%
  # 
  # # Take the string after the "p." to the end of string
  # mutate(VARIANT_P = str_match(VARIANT_c, ":p.(.*)$")[,2]) %>%
  # mutate(FUNCTION = str_match(INFO, "ExonicFunc.knownGene=(.*?);")[,2]) %>% # Or can do "ExonicFunc.ensGene="
  # mutate(COSMIC = str_match(INFO, "OccurSum=(.*?);")[,2]) %>%
  # mutate(ESP6500 = str_match(INFO, "esp6500siv2_all=(.*?);")[,2]) %>%

  # # separate FORMAT and DATA into format and data for corresponding value
  # mutate(format = FORMAT) %>%
  # mutate(read = DATA) %>%
  # # separate(format, into = paste("format", 1:10, sep="_"), # may not be necessary but can be useful
  # #          sep = ":") %>%
  # separate(read, into = paste("read", 1:10, sep="_"), # separate read aka DATA to be able to subset VAF and DEPTH
  #          sep = ":") %>%
  # mutate(VAF = read_3) %>% # genarate VAF and DEPTH var from read aka DATA
  # mutate(DEPTH = read_4) %>% # reorganize variable order
  select("patient_id", "CHROM", "POS", "REF", "ALT", "GENE", 
         "VARIANT_C", "VARIANT_P",
         "LOCATION", #"FUNCTION", 
         # "COSMIC", "ESP6500", "gnomAD_noncancer", "VAF",
         # "DEPTH", 
         "thousand_genome",
         "GT", "GQ","GQA","GTS","GTM","DP","AD","DP_Normal","AD_Normal",
         "INFO", "FORMAT", "DATA") %>%
  arrange(patient_id)

filtered_somatic_calls <- left_join(unique_patient_in_mutation, vcf6, by = "SLID_tumor")

saveRDS(filtered_somatic_calls, "filtered_somatic_calls.rds")

write_csv(filtered_somatic_calls, "filtered_somatic_calls.csv")


# Limit to the tumor SLID closest to blood
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "merging slid ID")
read_csv(paste0(path, "/List tumor SLID earliest or closest to germline.csv")) %>% 
  select(avatar_id, SLID_germline, collectiondt_germline, SLID_tumor, collectiondt_tumor) %>% 
  left_join(., filtered_somatic_calls, by = c("SLID_tumor", ))








