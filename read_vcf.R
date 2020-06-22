# Import library
library(tidyverse)
library(data.table)

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "read VCF files")
vcf <- read.delim(paste0(path, "/all_samples_sampleFiltered_bcfIsec.hg38_multianno_filtered_filtered_intersect.txt"),
                  na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))
# vcf <- read.delim(file.choose(), na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))

vcf <- vcf %>% 
  # Unite all mutations characteristics as a row names = make tidy data
  unite(col = "X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT",
        c("X.CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), sep = " ") %>% 
  # pivot longer which create 1 row per mutation per patients
  pivot_longer(-X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, 
               names_to = "patient_id", values_to = "DATA") %>% 
  # drop the mutations not present
  drop_na("DATA") %>% 
  # separate back the rownames to have mutation elements
  mutate(m_element = X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT) %>% 
  separate(m_element, 
           into = c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), 
           sep = " ") %>% 

# Generate the gene element variables from the INFO var
  mutate(GENE = str_match(INFO, "Gene.ensGene=(.*?);GeneDetail.ensGene=")[,2]) %>% 
  # VAVIANT_C
  # 1. Take the whole string before "esp6500siv2_all"
  mutate(VARIANT_C = str_match(INFO, "(.*);esp6500siv2_all")[,2]) %>% 
  mutate(VARIANT_c = VARIANT_C) %>% 
  # 2. Seaprate by "exon" and coalesce from last to first to get the last "exon" string
  separate(VARIANT_c, "exon", into = paste("VARIANT_c", 1:25, sep="_")) %>% 
  mutate(VARIANT_c = coalesce(VARIANT_c_25, VARIANT_c_24, VARIANT_c_23, VARIANT_c_22, VARIANT_c_21,
                              VARIANT_c_20, VARIANT_c_19, VARIANT_c_18, VARIANT_c_17, VARIANT_c_16,
                              VARIANT_c_15, VARIANT_c_14, VARIANT_c_13, VARIANT_c_12, VARIANT_c_11,
                              VARIANT_c_10, VARIANT_c_9, VARIANT_c_8, VARIANT_c_7, VARIANT_c_6,
                              VARIANT_c_5, VARIANT_c_4, VARIANT_c_3, VARIANT_c_2, VARIANT_c_1)) %>% 
  # 3. Take all between "c". and "p." but will lose data if "p." id not furnished
  mutate(VARIANT_C = str_match(VARIANT_c, "c.(.*):p.")[,2]) %>% 
  # 4. So do the same with just "c." and coalesce the ones lost (NA) by the "VARIANT_C1" which has them
  mutate(VARIANT_C1 = str_match(VARIANT_c, "c.(.*)")[,2]) %>% 
  mutate(VARIANT_C = coalesce(VARIANT_C, VARIANT_C1)) %>% 
  # VARIANT_P
  # Take the string after the "p." to the end of string
  mutate(VARIANT_P = str_match(VARIANT_c, ":p.(.*)$")[,2]) %>% 
  mutate(FUNCTION = str_match(INFO, "ExonicFunc.knownGene=(.*?);")[,2]) %>% # Or can do "ExonicFunc.ensGene="
  mutate(COSMIC = str_match(INFO, "OccurSum=(.*?);")[,2]) %>% 
  mutate(ESP6500 = str_match(INFO, "esp6500siv2_all=(.*?);")[,2]) %>%
  
  # separate FORMAT and DATA into format and data for corresponding value
  mutate(format = FORMAT) %>% 
  mutate(read = DATA) %>% 
  separate(format, into = paste("format", 1:10, sep="_"), # may not be necessary but can be useful
           sep = ":") %>% 
  separate(read, into = paste("read", 1:10, sep="_"), # separate read aka DATA to be able to subset VAF and DEPTH
           sep = ":") %>% 
    mutate(VAF = read_3) %>% # genarate VAF and DEPTH var from read aka DATA
    mutate(DEPTH = read_4) %>% # reorganize variable order
    select("patient_id", "CHROM", "POS", "REF", "ALT", "GENE", "VARIANT_C", 
           "VARIANT_P", "FUNCTION", "COSMIC", "ESP6500", "VAF", "DEPTH", "INFO", "FORMAT", "DATA") 


write_csv(vcf, paste0(path, "/cleaned vcf.csv"))

