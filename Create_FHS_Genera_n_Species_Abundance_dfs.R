"%&%" <- function(a,b) paste(a,b, sep = "")
library(tidyverse)
library(data.table)

#Use Phyloseq for diversity analysis
library(phyloseq)
#library(DECIPHER)
#library(phangorn)
library("microbiome")
library(ggplot2)
library(cowplot)
library(vegan)
"%ni%" <- Negate("%in%")

library(MMUPHin)
library(vegan)
library(magrittr)

library(microbiomeutilities)
library(ggdendro)
library(RColorBrewer)
library(metagenomeSeq)


#FHS

#Bring in the Meta Data
FHS_Meta <-read.table("/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Meta.txt")

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

#Keep only No-ANTIBIOTICS users and Covariates in the FHS_Meta
FHS_16s_rawcount <- subset(FHS_16s_rawcount, select = c(FHS_Meta$unique_ID))


FHStax <- as.matrix(FHStax)

#create Phyloseq Object with Tree
FHS_ps <- phyloseq(otu_table(FHS_16s_rawcount, taxa_are_rows=TRUE),
                   tax_table(FHStax), sample_data(FHS_Meta))
FHS_ps


#ALL FHS

fhs_glom_genus <- tax_glom(FHS_ps, taxrank = "Genus", NArm = T)
fhs_glom_species <- tax_glom(FHS_ps, taxrank = "Species", NArm = T)


#Convert to Genus relative abundance
fhs_glom <- microbiome::transform(fhs_glom_genus, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

rownames(fhs_totu) <- fhs_ttab$Genus
rownames(fhs_ttab) <- fhs_ttab$Genus

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Genus relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Genus.txt",
            quote = F, row.names = T, sep = "\t")
#Write out Genus Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Genus.txt",
            quote = F, row.names = T, sep = "\t")




#Convert to Species relative abundance
fhs_glom <- microbiome::transform(fhs_glom_species, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

#Create Microbes
Microbes <- vector(mode = "character", length = nrow(fhs_ttab))
# Genus_Specie
for (i in 1:length(Microbes)){
  Microbes[i] <- paste(fhs_ttab$Genus[i],fhs_ttab$Species[i], sep = "_")
}
fhs_ttab$Microbes <- Microbes


rownames(fhs_totu) <- fhs_ttab$Microbes
rownames(fhs_ttab) <- fhs_ttab$Microbes

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Genus relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Species.txt",
            quote = F, row.names = T, sep = "\t")

#Write out Genus Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Species.txt",
            quote = F, row.names = T, sep = "\t")






#MEN Only Abundance
men_meta <- subset(FHS_Meta, SEX==1)
men_abd <- subset(FHS_16s_rawcount, select = men_meta$unique_ID)

#Create Phyloseq Object
FHS_ps <- phyloseq(otu_table(men_abd, taxa_are_rows=TRUE),
                   tax_table(FHStax), sample_data(men_meta))
FHS_ps


fhs_glom_genus <- tax_glom(FHS_ps, taxrank = "Genus", NArm = T)
fhs_glom_species <- tax_glom(FHS_ps, taxrank = "Species", NArm = T)


#Convert Genus relative abundance
fhs_glom <- microbiome::transform(fhs_glom_genus, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

rownames(fhs_totu) <- fhs_ttab$Genus
rownames(fhs_ttab) <- fhs_ttab$Genus

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Genus relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Genus_Men_Only.txt",
            quote = F, row.names = T, sep = "\t")
#Write out Genus Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Genus_Men_Only.txt",
            quote = F, row.names = T, sep = "\t")




#Convert Species relative abundance
fhs_glom <- microbiome::transform(fhs_glom_species, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

#Create Microbes
Microbes <- vector(mode = "character", length = nrow(fhs_ttab))
# Genus_Specie
for (i in 1:length(Microbes)){
  Microbes[i] <- paste(fhs_ttab$Genus[i],fhs_ttab$Species[i], sep = "_")
}
fhs_ttab$Microbes <- Microbes


rownames(fhs_totu) <- fhs_ttab$Microbes
rownames(fhs_ttab) <- fhs_ttab$Microbes

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Species relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Species_Men_Only.txt",
            quote = F, row.names = T, sep = "\t")

#Write out Species Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Species_Men_Only.txt",
            quote = F, row.names = T, sep = "\t")










#WOMEN Only Abundance
women_meta <- subset(FHS_Meta, SEX==2)
women_abd <- subset(FHS_16s_rawcount, select = women_meta$unique_ID)

#Create Phyloseq Object
FHS_ps <- phyloseq(otu_table(women_abd, taxa_are_rows=TRUE),
                   tax_table(FHStax), sample_data(women_meta))
FHS_ps


fhs_glom_genus <- tax_glom(FHS_ps, taxrank = "Genus", NArm = T)
fhs_glom_species <- tax_glom(FHS_ps, taxrank = "Species", NArm = T)


#Convert Genus relative abundance
fhs_glom <- microbiome::transform(fhs_glom_genus, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

rownames(fhs_totu) <- fhs_ttab$Genus
rownames(fhs_ttab) <- fhs_ttab$Genus

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Genus relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Genus_Women_Only.txt",
            quote = F, row.names = T, sep = "\t")
#Write out Genus Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Genus_Women_Only.txt",
            quote = F, row.names = T, sep = "\t")




#Convert Species relative abundance
fhs_glom <- microbiome::transform(fhs_glom_species, "compositional")

fhs_ttab <- tax_table(fhs_glom) %>% as.data.frame()
fhs_totu <- otu_table(fhs_glom) %>% as.data.frame()

#Create Microbes
Microbes <- vector(mode = "character", length = nrow(fhs_ttab))
# Genus_Specie
for (i in 1:length(Microbes)){
  Microbes[i] <- paste(fhs_ttab$Genus[i],fhs_ttab$Species[i], sep = "_")
}
fhs_ttab$Microbes <- Microbes


rownames(fhs_totu) <- fhs_ttab$Microbes
rownames(fhs_ttab) <- fhs_ttab$Microbes

#Transpose OTU so that column is microbes, and rows are samples
fhs_totu <- t(fhs_totu) %>% as.data.frame()

#Write out Species relative abundance OTU and Taxa
write.table(fhs_totu,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Abd_Species_Women_Only.txt",
            quote = F, row.names = T, sep = "\t")

#Write out Species Taxa
write.table(fhs_ttab,file = "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/O2_FHS_16S_Taxa_Species_Women_Only.txt",
            quote = F, row.names = T, sep = "\t")




