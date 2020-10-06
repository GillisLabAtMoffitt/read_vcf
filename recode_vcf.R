# Import library
library(tidyverse)
library(data.table)

# Load data
vcf <- read.delim(file.choose()) %>% 
  mutate_all(funs(str_replace(., "\\.\\/(.*)", NA_character_)))


# Data cleaning
vcf <- vcf %>% 
  # Unite all mutations characteristics as a row names = make tidy data to allow pivot
  unite(col = "X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT",
        c("X.CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), sep = " ") %>% 
  # pivot longer which create 1 row per mutation per patients
  pivot_longer(-X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, 
               names_to = "patient_id", values_to = "DATA") 

# Save each patient ID with or without mutation to add back at the end
unique_patient_in_mutation <- as.data.frame(unique(vcf$patient_id)) %>% 
  `colnames<-` (c("patient_id"))

vcf <- vcf %>% 
  # drop the mutations not present (make code faster)
  drop_na("DATA") %>% 
  # separate back the rownames to have mutation elements
  mutate(m_element = X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT) %>% 
  separate(m_element, 
           into = c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), 
           sep = " ") %>% 
  # Generate the gene element variables from the INFO var
  mutate(GENE = str_match(INFO, ";Gene.knownGene=(.*?);")[,2]) %>% # was before "Gene.ensGene=(.*?);GeneDetail.ensGene="
  
  # VAVIANT_C
  # 1. Take the whole string before "esp6500siv2_all"
  mutate(VARIANT_C = str_match(INFO, "(.*);esp6500siv2_all")[,2]) %>% 
  mutate(VARIANT_c = VARIANT_C) %>% 
  # 2. Seaprate by "exon" and coalesce from last to first to get the last "exon" string
  mutate(VARIANT_c = str_split(string= VARIANT_c, pattern= "exon")) %>%
  unnest(VARIANT_c) %>%
  group_by(patient_id, VARIANT_C) %>%
  mutate(row = row_number()) %>%
  spread(row, VARIANT_c) %>% 
  ungroup() %>% 
  # separate(VARIANT_c, "exon", into = paste("VARIANT_c", 1:25, sep="_")) %>% 
  mutate(VARIANT_c = coalesce(!!! select(., last_col():"1"))) %>%
  # 3. Take all between "c". and "p." but will lose data if "p." id not furnished
  mutate(VARIANT_C = str_match(VARIANT_c, "c.(.*):p.")[,2]) %>% 
  # 4. So do the same with just "c." and coalesce the ones lost (NA) by the "VARIANT_C1" which has them
  mutate(VARIANT_C1 = str_match(VARIANT_c, "c.(.*)")[,2]) %>% 
  mutate(VARIANT_C = coalesce(VARIANT_C, VARIANT_C1)) %>% 
  
  # VARIANT_P
  # Take the string after the "p." to the end of string
  mutate(VARIANT_P = str_match(VARIANT_c, ":p.(.*)$")[,2]) %>% 
  mutate(LOCATION = str_match(INFO, ";Func.ensGene=(.*?);")[,2]) %>% 
  mutate(FUNCTION = str_match(INFO, "ExonicFunc.knownGene=(.*?);")[,2]) %>% # Or can do "ExonicFunc.ensGene="
  mutate(COSMIC = str_match(INFO, "OccurSum=(.*?);")[,2]) %>% 
  mutate(ESP6500 = str_match(INFO, "esp6500siv2_all=(.*?);")[,2]) %>%
  
  # separate FORMAT and DATA into format and data for corresponding value
  mutate(format = FORMAT) %>% 
  mutate(read = DATA) %>% 
  # separate(format, into = paste("format", 1:10, sep="_"), # may not be necessary but can be useful
  #          sep = ":") %>% 
  separate(read, into = paste("read", 1:10, sep="_"), # separate read aka DATA to be able to subset VAF and DEPTH
           sep = ":") %>% 
    mutate(VAF = read_3) %>% # genarate VAF and DEPTH var from read aka DATA
    mutate(DEPTH = read_4) %>% # reorganize variable order
    select("patient_id", "CHROM", "POS", "REF", "ALT", "GENE", "VARIANT_C", 
           "VARIANT_P", "LOCATION", "FUNCTION", "COSMIC", "ESP6500", "VAF", 
           "DEPTH", "INFO", "FORMAT", "DATA") %>% 
  arrange(patient_id)



vcf <- left_join(unique_patient_in_mutation, vcf, by = "patient_id")

# Saving data
# To save in the working directory
write_csv(vcf, "cleaned vcf.csv")
# To save in an other working directory
# Set up the directory first, here my desktop
setwd("/Users/colinccm/Desktop")
write_csv(vcf, "cleaned vcf.csv")
