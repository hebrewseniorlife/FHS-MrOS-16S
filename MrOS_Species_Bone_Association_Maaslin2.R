"%&%" <- function(a,b) paste(a,b, sep = "")
library(tidyverse)
library(MMUPHin)
library(vegan)
library(magrittr)
library(Maaslin2)
"%ni%" <- Negate("%in%")

#"S:/Obesity_Bone/Paul/Data/data_to_O2/"

data_path <- "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/"
res_path <- "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/16S_results/SPECIES/MrOS/"

#Load data created on #HSL RStudio
meta <- read.table(file = data_path %&% "O2_MrOS_16S_Meta.txt", 
                   sep = "\t",header = T, row.names = 1, check.names = F)


#Create MrOS Clinic Site
clinic <- vector(mode="character", length = nrow(meta))
for (i in 1:nrow(meta)) {
  x1 <- strsplit(meta$ID[i], "", fixed = T)[[1]][1]
  x2 <- strsplit(meta$ID[i], "", fixed = T)[[1]][2]
  clinic[i] <- paste(x1,x2,sep = "")
}
meta$CLINIC <- clinic
meta$CLINIC <- as.factor(meta$CLINIC)
# table(clinic)
# clinic
# BI  MN  PA  PI  PO  SD 
# 104 144 109 122 156 201 



#Bring
#Bring in the Batches
batches <-read.table(data_path %&% "MrOS_The_Missing_320_Samples.txt",
                     header = T)

#RENAME MrOS SAMPLE ID
x <- batches$ID
#Split and bring out the MROS UNIQUE ID similar to ID in cov and hrpqct
for (i in 1:length(x)){
  x[i] <- strsplit(x[i], split = "_", fixed = T)[[1]][2]
}
batches$ID <- x

#Split metadata into Batches
b1 <- subset(meta, ID %in% batches$ID)
b1$group <- "Batch1"
b2 <- subset(meta, ID %ni% batches$ID)
b2$group <- "Batch2"

b1b2 <- rbind(b1, b2)
meta <- b1b2
meta$group <- as.factor(meta$group)


#Do Batch Adjustment with MMUPHIN

#######################       START RUNNING MMUPHIN !       #################
#https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html

covs <- c("SEX", "AGE", "SMOKING", "BMI", "DIABETES", "HOSPITALIZATION", "TOTAL_MEDS",
          "PPI", "METFORMIN", "RACE", "DQI")

x <- t(abd) %>% data.frame() #make sample names to be columns
#Keep sample names in the same order
x <- subset(x, select = c(meta$ID))

fit_adjust_batch <- adjust_batch(feature_abd = x,
                                 batch = "group",
                                 covariates = NULL,
                                 data = meta,
                                 control = list(verbose = FALSE,
                                                output=res_path))

adj_abd <- fit_adjust_batch$feature_abd_adj

#Check the effect of the batch adjustment to assess the total variability in microbial
#profiles attributable to batch differences using PERMANOVA

D_before <- vegdist(t(x))
D_after <- vegdist(t(adj_abd))

set.seed(123)
fit_adonis_before <- adonis(D_before ~ group,
                            data = meta, permutations=999, method = "bray")
fit_adonis_after <- adonis(D_after ~ group,
                           data = meta, permutations=999, method = "bray")
#Before
print(fit_adonis_before)
#After
print(fit_adonis_after)
df <- fit_adonis_after$aov.tab
df

#Remove batch from the metadata df
meta$group <- NULL


#Use Adjusted Abundance table for Subsequent analysis
#Make microbes to be the columns
adj_abd <- t(adj_abd) %>% data.frame()

#Write out adjusted abundance
write.table(adj_abd, data_path %&% "adj_O2_MrOS_16S_Abd_Species.txt",
            quote = F, sep = "\t", row.names = T)





#Load adjusted abundance data created on #O2 RStudio, with this script above
adj_abd <- read.table(file = data_path %&% "adj_O2_MrOS_16S_Abd_Species.txt", 
                      sep = "\t", header = T, row.names = 1, check.names = F)



#RUN MAASLIN2 https://forum.biobakery.org/t/maaslin2-heatmap-vs-significance/1800/3
#https://groups.google.com/g/picrust-users/c/XwaHlrittOQ/m/A2PQz_XSAgAJ

bone_pheno <- c("t_Ct_Th","t_Ct_vBMD","t_CtBATA","t_Tb_Inn_vBMD","t_Tb_N",
                "t_Tb_vBMD","t_Tot_vBMD","t_Total_Area","t_fail_load_FEA","r_Ct_Th",
                "r_Ct_vBMD","r_CtBATA","r_Tb_Inn_vBMD","r_Tb_N","r_Tb_vBMD",
                "r_Tot_vBMD","r_Total_Area","r_fail_load_FEA")

covs <- c("AGE", "SMOKING", "BMI", "DIABETES", "HOSPITALIZATION", "TOTAL_MEDS",
          "PPI", "METFORMIN", "RACE", "DQI")

#USE STANDARDIZED HRPQCT
#Mean=0, SD=1
x_meta <- meta
x <- subset(x_meta, select = bone_pheno)
x <- x %>% mutate_at(bone_pheno, ~(scale(.) %>% as.vector))
#remove the bone_pheno, and keep the covariates
x_meta[,bone_pheno] <- NULL
#join the standardized hrpqct with covariates. they are already in same sample order
x_meta <- cbind(x, x_meta)



#LM
#RUN MAASLIN2 with this settings. That is "LM" and "TSS", transform LOG

for (i in 1:length(bone_pheno)){
  x <- c(bone_pheno[i],covs) #create the list of fixed effects

  fit_fhs_maaslin2 <- Maaslin2(input_data = adj_abd, input_metadata = x_meta,
                               output = res_path %&% "HRPQCT/" %&% bone_pheno[i],
                               fixed_effects = x, analysis_method = "LM",
                               normalization = "TSS", transform = "LOG",
                               max_significance = 0.25,
                               random_effects = "CLINIC")
}



