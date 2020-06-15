# Import library
library(tidyverse)
library(data.table)

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "read VCF files")
vcf <- read.delim(paste0(path, "/all_samples_sampleFiltered_bcfIsec.hg38_multianno_filtered_filtered_intersect.txt"),
                  na.strings = c("./.:.:.:.:.:.:.", "./.:.:.:.:.:.:.:.:.:."))


