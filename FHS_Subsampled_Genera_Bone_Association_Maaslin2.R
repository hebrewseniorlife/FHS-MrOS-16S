"%&%" <- function(a,b) paste(a,b, sep = "")
library(tidyverse)
library(Maaslin2)

#data_path <- "S:/Obesity_Bone/Paul/Data/data_to_O2/"

data_path <- "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/data/"
res_path <- "/n/groups/ifaromics/paul/HRPQCT_16S_Analysis/16S_results/GENUS/subsample_FHS/corrected_seed_1234/"

#Load data created on #HSL RStudio
abd <- read.table(file = data_path %&% "O2_FHS_Subsampled_16S_Abd_Genus.txt", 
                  sep = "\t", header = T, row.names = 1, check.names = F)
meta <- read.table(file = data_path %&% "O2_FHS_16S_Meta.txt", 
                   sep = "\t",header = T, row.names = 1, check.names = F)

#Subset metadata to contain only the sub-sampled samples
meta <- subset(meta, rownames(meta) %in% rownames(abd))

#RUN MAASLIN2 https://forum.biobakery.org/t/maaslin2-heatmap-vs-significance/1800/3
#https://groups.google.com/g/picrust-users/c/XwaHlrittOQ/m/A2PQz_XSAgAJ

bone_pheno <- c("t_Ct_Th","t_Ct_vBMD","t_CtBATA","t_Tb_Inn_vBMD","t_Tb_N",
                "t_Tb_vBMD","t_Tot_vBMD","t_Total_Area","t_fail_load_FEA","r_Ct_Th",
                "r_Ct_vBMD","r_CtBATA","r_Tb_Inn_vBMD","r_Tb_N","r_Tb_vBMD",
                "r_Tot_vBMD","r_Total_Area","r_fail_load_FEA")

covs <- c("SEX", "AGE", "SMOKING", "BMI", "DIABETES", "HOSPITALIZATION", "TOTAL_MEDS",
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
  
  fit_fhs_maaslin2 <- Maaslin2(input_data = abd, input_metadata = x_meta,
                               output = res_path %&% "HRPQCT/" %&% bone_pheno[i],
                               fixed_effects = x, analysis_method = "LM",
                               normalization = "TSS", transform = "LOG",
                               max_significance = 0.25)
}


