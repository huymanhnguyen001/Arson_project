# Loading Packages --------------------------------------------------------
library(ggplot2)
library(readxl)
library(ggpubr)
library(gtable)
library(gridExtra)
library(viridis)
library(wesanderson)
library(tidyverse)
library(lubridate)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(data.table)
library(corrplot)

# Analysis Progress -------------------------------------------------------
# fail to appropriately pivot_wider: each column should be name of compound and each row will be
# the measurement of that compount, such as MF, RMF, RT1, RT2, Area %, Height, Ion 1, Ion 2, etc.

# Excel file name explanation ---------------------------------------------------
# The attached files are direct injection IL and the following explains the file naming:
# - 0220F00991
# - month (02) year (20) station (009) octane (91)
# Anything with a D at end is diesel and the composites are obviously a mixture of a fuel. 


# Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')
df_list <- map(file_list, read_excel, sheet = 'Results')

# Combine all data
all_data <- bind_rows(df_list)

# Filtering out column bleed and solvent --------------------------------------
filter_list <- c("Carbon disulfide",
                 "Benzene",
                 "Cyclotrisiloxane",
                 "Toluene",
                 "Ethylbenzene",
                 "Xylene")

for (ele in filter_list) {
  all_data <- all_data %>%
    filter(!grepl(ele, Compound))
}

# Number of unique compounds after filering
length(unique(all_data$Compound)) # 9505

# Checkpoint for dimethylbenzene
str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene

# Distribution of Peak Area and Peak Height values
hist(all_data$Area, xlim = c(0, 24945800)) # 95% of peak area is 0 - 2E+07
hist(all_data$Height, xlim = c(0, 389655)) # 95% of peak height is 0 - 2E+05

# Scaling Peak Area and Peak Height
all_data$Area <- scale(all_data$Area)
all_data$Height <- scale(all_data$Height)

# Shapiro-Wilk Normality Test ---------------------------------------------
# Peak Area

# Peak Height


# Data Normalization -----------------------------------------------------


# summarize similar compound detection by the mean value  
for (i in 1:length(df_list)) {
  df_list[[i]] <- df_list[[i]] %>% 
    group_by(Compound) %>% 
    summarise(across(everything(), mean))#  %>%
    # mutate(sample_name = file_list[[i]]) %>%
    # relocate(sample_name) %>%
    # mutate(Compound_and_sample = paste(Compound, "-", sample_name)) %>%
    # relocate(Compound_and_sample)
}


# Convert list of excel file (tibble) into a single data (tibble) ---------
all_data <- rbindlist(df_list)


# Pivot wider for PCA -----------------------------------------------------
df_list[[1]] %>% 
  pivot_longer(-c(Compound_and_sample, sample_name, Compound)) %>%
  pivot_wider(names_from = Compound, values_from = value)

# i <- 1
# df_list_wider <- list()
# df_wider$sample_name <- file_list[[i]]
# df_wider <- relocate(df_wider, sample_name)
# df_list_wider[[i]] <- df_wider
# i <- i + 1


# PCA ---------------------------------------------------------------------
file1 <- copy(df_list[[1]])
file1 <- file1 %>%
  column_to_rownames(., var = "Compound")
file1_pca <- PCA(file1, scale.unit = TRUE, ncp = 5, graph = FALSE)

fviz_eig(file1_pca,
         addlabels = TRUE) 

var <- get_pca_var(file1_pca)
corrplot(var$cos2, is.corr=FALSE)
fviz_cos2(file1_pca, choice = "var", axes = 1:2)

ind <- get_pca_ind(file1_pca)
fviz_cos2(file1_pca, choice = "ind")
fviz_contrib(file1_pca, choice = "ind", axes = 1:2)
