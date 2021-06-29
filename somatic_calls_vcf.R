# Import library
library(tidyverse)
library(data.table)
library(gtsummary)

doParallel::registerDoParallel()
# Load data
# data <- read_delim(file.choose(), delim = "\t")
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
unique_patient_in_mutation <- as.data.frame(unique(vcf$SLID_tumor)) %>% 
  `colnames<-` (c("SLID_tumor"))


vcf1 <- vcf %>% 
  # drop the mutations not present (make code faster)
  drop_na("DATA")

saveRDS(vcf1, "mid processed data.rds")

vcf2 <- vcf1 %>% 
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
  
  select("SLID_tumor", "CHROM", "POS", "REF", "ALT", "GENE", 
         "VARIANT_C", "VARIANT_P",
         "LOCATION", #"FUNCTION", 
         # "COSMIC", "ESP6500", "gnomAD_noncancer", "VAF",
         # "DEPTH", 
         "thousand_genome",
         "GT", "GQ","GQA","GTS","GTM","DP","AD","DP_Normal","AD_Normal",
         "INFO", "FORMAT", "DATA") %>%
  arrange(SLID_tumor)

filtered_somatic_calls <- left_join(unique_patient_in_mutation, vcf6, by = "SLID_tumor")

saveRDS(filtered_somatic_calls, "filtered_somatic_calls.rds")

write_csv(filtered_somatic_calls, "filtered_somatic_calls.csv")


# Limit to the tumor SLID closest to blood
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "merging slid ID")
final <- read_csv(paste0(path, "/List tumor SLID earliest or closest to germline.csv")) %>% 
  select(avatar_id, SLID_germline, collectiondt_germline, SLID_tumor, collectiondt_tumor) %>% 
  left_join(., filtered_somatic_calls, by = "SLID_tumor")

write_csv(final, "final.csv")


################ MAF file

maf_file <- read_delim(file.choose(), delim = "\t")

final <- read_csv(paste0(path, "/List tumor SLID earliest or closest to germline.csv")) %>% 
  select(avatar_id, SLID_germline, collectiondt_germline, SLID_tumor, collectiondt_tumor) %>% 
  left_join(., maf_file, by = c("SLID_tumor" = "Tumor_Sample_Barcode"))

final <- final %>% 
  mutate(tumor_VAF = t_alt_count / t_depth) %>% 
  mutate(normal_VAF = n_alt_count / n_depth) %>% 
  purrr::keep(~!all(is.na(.))) %>% 
  
  select("avatar_id", "SLID_germline", "collectiondt_germline", "SLID_tumor", 
         "collectiondt_tumor", "Hugo_Symbol", "NCBI_Build", "Chromosome",
         "Start_Position", "End_Position", "Strand", 
         "Variant_Classification", "Variant_Type", "Reference_Allele", 
         "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
         "HGVSc", "HGVSp_Short", "Transcript_ID", "Exon_Number", 
         "tumor_VAF", "normal_VAF",
         "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count",
         "n_alt_count", everything())

write_csv(final, "somatic mutation.csv")

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", 
                 "CHIP in Avatar", "TumorMuts")

tbl <- final %>% 
  distinct(avatar_id, Hugo_Symbol, .keep_all = TRUE) %>% 
  select(Hugo_Symbol) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  as_gt() %>% 
  gt::tab_source_note(gt::md("**Count 1 gene per patient**"))

gt::gtsave(tbl, zoom = 1, paste0(path, "/Output/Tumor Mutations in MM Avatar.pdf"))

tbl <- final %>% 
  distinct(avatar_id, Hugo_Symbol, .keep_all = TRUE) %>% 
  filter(str_detect(Hugo_Symbol, "KRAS|NRASTRAF3|DIS3|BRAF|FAM46C|TP53|CYLD|MAX|RB1")) %>%
  select(Hugo_Symbol) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  as_gt() %>% 
  gt::tab_source_note(gt::md("**Count 1 gene per patient**"))

gt::gtsave(tbl, zoom = 1, paste0(path, "/Output/Tumor Mutations MM genes.pdf"))

list <- final %>% 
  select(Variant_Classification) %>% count(Variant_Classification, sort = TRUE)

write_csv(list, paste0(path, "/Output/List funtions in data.csv"))  
  
list <- final %>% 
  filter(Hugo_Symbol != "Unknown") %>%
  select(Variant_Classification) %>% 
  count(Variant_Classification, sort = TRUE)

write_csv(list, paste0(path, "/Output/List funtions in data Unknown gene removed.csv")) 
  
list <- final %>% 
  select(Exon_Number) %>% 
  count(Exon_Number, sort = TRUE)

write_csv(list, paste0(path, "/Output/List exon number in data.csv")) 
  
  
  
  
  
  
  

