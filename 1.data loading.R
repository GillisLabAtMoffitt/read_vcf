# Import library
library(tidyverse)
library(data.table)

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "read VCF files")
vcf <- read.delim(paste0(path, "/all_samples_sampleFiltered_bcfIsec.hg38_multianno_filtered_filtered_intersect.txt"),
                  na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))

# Unite all mutations characteristics as a row names
colnames(vcf)
vcf1 <- vcf %>%
  unite(col = "X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT",
        c("X.CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), sep = " * ")

# All us to pivot longer which create 1 row per mutation per patients, plus drop the mutations not found
vcf2 <- vcf1 %>%
  pivot_longer(-X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, 
               names_to = "patient_id", values_to = "DATA") %>% 
  drop_na("DATA")

vcf3 <- vcf2 %>% 
  mutate(m_element = X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT) %>% 
  separate(m_element, 
           into = c("X.CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), 
           sep = " ")

vcf3$GENE <- str_match(vcf3$INFO, "Gene.ensGene=(.*?);GeneDetail.ensGene=")[,2]
vcf3$VARIANT_C <- str_match(vcf3$INFO,"c.(.*?);esp6500siv2_all")
vcf3$VARIANT_P <- str_match(vcf3$INFO, "p.(.*?);esp6500siv2_all")
vcf3$FUNCTION <- str_match(vcf3$INFO, "ExonicFunc.knownGene=(.*?);")[,2] # Or can do "ExonicFunc.ensGene="
vcf3$COSMIC <- str_match(vcf3$INFO, "OccurSum=(.*?);")[,2]
vcf3$ESP6500 <- str_match(vcf3$INFO, "esp6500siv2_all=(.*?);")[,2]

vcf4 <- vcf3 %>% 
  select("patient_id", "X.CHROM", "POS",  "REF", "ALT", "INFO", "FORMAT", "DATA") %>% # remove "ID", "QUAL"
  mutate(read)
  separate(FORMAT, into = paste("FORMAT", 1:10, sep="_"), 
           sep = ":") %>% 
  separate(DATA, into = paste("read", 1:10, sep="_"),
           sep = ":") 

