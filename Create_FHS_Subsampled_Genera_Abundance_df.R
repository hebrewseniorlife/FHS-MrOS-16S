"%&%" <- function(a,b) paste(a,b, sep = "")
library(tidyverse)
library(data.table)
library(phyloseq)
library("microbiome")
library(ggplot2)
library(cowplot)
library(vegan)
"%ni%" <- Negate("%in%")

library(microbiomeutilities)
library(ggdendro)
library(RColorBrewer)
library(metagenomeSeq)

#FHS

#Bring in the Meta Data
FHS_Meta <-read.table("/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Meta.txt")
#FHS_Meta$Cohort <- "FHS"

#Bring in the 16S Sequence data
FHS_16s_rawcount <- read.table("/n/groups/ifaromics/paul/picrust2_data/FHS_workflows_dada2_silva_maxee2_maxmismatch3/all_samples_taxonomy_closed_reference_taxcolumns.tsv", header=T,
                               sep = "\t", check.names = F, row.names = 1)
#takeout taxonomy
FHStax <- FHS_16s_rawcount[,c("Kingdom", "Phylum", "Class", "Order", "Family", 
                              "Genus", "Species")]
#remove taxonomy
FHS_16s_rawcount[,c("Kingdom", "Phylum", "Class", "Order", "Family", 
                    "Genus", "Species")] <- NULL

#Note, some of the ID_shorts were changed during the dada2 analysis such that special
#character underscore (_) were changed to dash "-"
#so change it back here
fhs_16s_col <- data.frame(old_ID_short=colnames(FHS_16s_rawcount))
fhs_16s_col$correct_ID <- ""
for (i in 1:nrow(fhs_16s_col)) { #the change will happen only where there's match. no error
  fhs_16s_col$correct_ID[i] <- sub("-","_", fhs_16s_col$old_ID_short[i], fixed = TRUE)
}

#retain only samples with no bowel surgery or use of antibiotics
#FHS
fhsnobiotics <-  read.table("/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/fhs_dada2_16s_nobowelsurgery_antibiotic_use.txt",
                            header = T, sep = "\t")
fhs_16s_col <- inner_join(fhs_16s_col, fhsnobiotics, by = c("correct_ID"="ID_short"))

#change colnames of the 16s data from old_ID_short to the FHS unique_ID
FHS_16s_rawcount <- subset(FHS_16s_rawcount, select = c(fhs_16s_col$old_ID_short))
colnames(FHS_16s_rawcount) <- c(fhs_16s_col$unique_ID)


#Get the exact samples in the previous FHS Abd Table
data_path <- "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/"
#Load data created on #HSL RStudio
abd <- read.table(file = data_path %&% "O2_FHS_16S_Abd_Genus.txt", 
                  sep = "\t", header = T, row.names = 1, check.names = F)
#SUBSAMPLE out 836 samples, randomly
fhs_samples <- rownames(abd)
#Set Seed for replicability
set.seed(1234)
fhs_800 <- fhs_samples[sample(1:length(fhs_samples), size = 836, replace = F)]



#Subset abd and meta with subsamples
#abd <- subset(abd, rownames(abd) %in% fhs_800)
FHS_Meta <- subset(FHS_Meta, rownames(FHS_Meta) %in% fhs_800)
table(FHS_Meta$SEX)
#Keep only the subsampled samples in the FHS_Meta
FHS_16s_rawcount <- subset(FHS_16s_rawcount, select = c(FHS_Meta$unique_ID))


FHStax <- as.matrix(FHStax)

#create Phyloseq Object with Tree
FHS_ps <- phyloseq(otu_table(FHS_16s_rawcount, taxa_are_rows=TRUE),
                   tax_table(FHStax), sample_data(FHS_Meta))
FHS_ps

#Glom at Genus
fhs_glom_genus <- tax_glom(FHS_ps, taxrank = "Genus", NArm = T)


#Convert to Genus relative abundance
fhs_glom <- microbiome::transform(fhs_glom_genus, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

rownames(fhs_totu) <- fhs_ttab$Genus
rownames(fhs_ttab) <- fhs_ttab$Genus

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Genus relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_Subsampled_16S_Abd_Genus.txt",
            quote = F, row.names = T, sep = "\t")
#Write out Genus Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_Subsampled_16S_Taxa_Genus.txt",
            quote = F, row.names = T, sep = "\t")




