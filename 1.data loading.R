# Import library
library(tidyverse)
library(data.table)

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "read VCF files")
vcf <- read.delim(paste0(path, "/all_samples_sampleFiltered_bcfIsec.hg38_multianno_filtered_filtered_intersect.txt"),
                  na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))
vcf <- read.delim(file.choose(), na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))

vcf <- vcf %>% 
  # Unite all mutations characteristics as a row names
  unite(col = "X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT",
        c("X.CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), sep = " ") %>% 
  # pivot longer which create 1 row per mutation per patients
  pivot_longer(-X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT, 
               names_to = "patient_id", values_to = "DATA") %>% 
  # drop the mutations not found
  drop_na("DATA") %>% 
  # re-separate the rownames to have mutation elements
  mutate(m_element = X.CHROM_POS_ID_REF_ALT_QUAL_FILTER_INFO_FORMAT) %>% 
  separate(m_element, 
           into = c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), 
           sep = " ")
# Generate the gene element variables fron the INFO var
vcf$GENE <- str_match(vcf$INFO, "Gene.ensGene=(.*?);GeneDetail.ensGene=")[,2]
vcf$VARIANT_C <- str_match(vcf$INFO, "(.*);esp6500siv2_all")[,2]
vcf$VARIANT_C
vcf$VARIANT_C <- str_sub(vcf$VARIANT_C, start = -60)
vcf$VARIANT_C
vcf$VARIANT_C <- str_match(vcf$VARIANT_C, "c.(.*)")[,2]
vcf$VARIANT_C
vcf$VARIANT_P <- str_match(vcf$VARIANT_C, ":p.(.*)$")
vcf$VARIANT_P
vcf$VARIANT_C <- str_match(vcf$VARIANT_C, "(.*):p.")[,1]
vcf$VARIANT_C
vcf$FUNCTION <- str_match(vcf$INFO, "ExonicFunc.knownGene=(.*?);")[,2] # Or can do "ExonicFunc.ensGene="
vcf$COSMIC <- str_match(vcf$INFO, "OccurSum=(.*?);")[,2]
vcf$ESP6500 <- str_match(vcf$INFO, "esp6500siv2_all=(.*?);")[,2]

vcf <- vcf %>% # separate FORMAT and DATA into format and data for corresponding value
  mutate(format = FORMAT) %>% 
  mutate(read = DATA) %>% 
  separate(format, into = paste("format", 1:10, sep="_"), # may not be necessary but can be useful
           sep = ":") %>% 
  separate(read, into = paste("read", 1:10, sep="_"), # separate read aka DATA to be able to subset VAF and DEPTH
           sep = ":") %>% 
    mutate(VAF = read_3) %>% # genarate VAF and DEPTH var from read aka DATA
    mutate(DEPTH = read_4) %>% # reorganize variable order
    select("patient_id", "CHROM", "POS", "REF", "ALT", "GENE", "VARIANT_C", 
           "VARIANT_P", "FUNCTION", "COSMIC", "ESP6500", "VAF", "DEPTH", "INFO", "FORMAT", "DATA",
           "format_1", "read_1", "format_2", "read_2", "format_3", "read_3", "format_4", "read_4", 
           "format_5", "read_5", "format_6", "read_6", "format_7", "read_7", "format_8", "read_8", 
           "format_9", "read_9", "format_10", "read_10") 


write_csv(vcf, paste0(path, "/vcf.csv"))

