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

set.seed(12345)
# Analysis Progress -------------------------------------------------------
# fail to appropriately pivot_wider: each column should be name of compound and each row will be
# the measurement of that compount, such as MF, RMF, RT1, RT2, Area %, Height, Ion 1, Ion 2, etc.

# Excel file name explanation ---------------------------------------------------
# The attached files are direct injection IL and the following explains the file naming:
# - 0220F00991
# - month (02) year (20) station (009) octane (91)
# Anything with a D at end is diesel and the composites are obviously a mixture of a fuel. 


# Functions -------------------------------------------------------------------------------------------------------
# Filtering matched compound names
filtering <- function (data, filter_list) { #, percentile, column_list
  clean_data <- copy(data)
  for (ele in filter_list) {
    clean_data <- clean_data %>%
      filter(!grepl(ele, Compound))
  }
  return(clean_data)
}

# Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')
df_list <- map(file_list, read_excel, sheet = 'Results')

# Labeling duplicate compound names with numeric suffix and with sample name 
# for (i in 1:length(df_list)) {
#   df_list[[i]] <- df_list[[i]] %>% 
#     group_by(Compound) %>%
#     mutate(Compound = make.names(Compound, unique = TRUE, allow_ = FALSE)) %>%
#     mutate(Compound = paste(Compound, "-", file_list[[i]]))
# }

# Combine all data
# all_data <- bind_rows(df_list)

# Checkpoint for dimethylbenzene
# str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene


# Filtering out column bleed and solvent --------------------------------------
# \< & \>: "\<" is an escape sequence for the beginning of a word, and ">" is used for end
# "\b" is an anchor to identify word before/after pattern

df_list_clean <- map(df_list, filtering, filter_list = c("Carbon.disulfide",
                                        "Benzene - ", 
                                        "Cyclotrisiloxane..hexamethyl",
                                        "Cyclotetrasiloxane..octamethyl",
                                        "Toluene",
                                        "Ethylbenzene",
                                        "Xylene")) 

df_list_clean_1 <- df_list_clean[[1]]

df_list_clean_1 <- df_list_clean_1 %>%
  select(-ends_with(c("Area %", "Ion 1", "Ion 2", "Ion 3")))
# Calculate percentage Peak Area and Peak Height
df_list_clean_1$Percent_Area <- df_list_clean_1$Area/sum(df_list_clean_1$Area)
df_list_clean_1$Percent_Height <- df_list_clean_1$Height/sum(df_list_clean_1$Height)

# Plotting the Percentage area/height vs. compound
df_list_clean_1 <- df_list_clean_1 %>%
  arrange(Percent_Height, Percent_Area)


# Filter peak area based on percentile  ---------------------------------------------------------------------------
# check percentile distribution
quantile(all_data_clean$Area) #summary()
quantile(all_data_clean$Height) #summary()

# filter iteration --> min. cut = mean(min., 25th percentile), max. cut = mean(max., 75th percentile)
# while (nrow(all_data_clean) > 4000) {
#   filter_quantile <- subset(all_data_clean, (Area > mean(quantile(all_data_clean$Area)[1],
#                                                          quantile(all_data_clean$Area)[2])
#                                              & Area < mean(quantile(all_data_clean$Area)[4],
#                                                            quantile(all_data_clean$Area)[5])) &
#                               (Height > mean(quantile(all_data_clean$Height)[1],
#                                              quantile(all_data_clean$Height)[2])
#                                & Height < mean(quantile(all_data_clean$Height)[4],
#                                                quantile(all_data_clean$Height)[5])))
#   new_area_quantile <- quantile(filter_quantile$Area)
#   new_height_quantile <- quantile(filter_quantile$Height)
#   new_filter_quantile <- subset(all_data_clean, (Area > mean(quantile(all_data_clean$Area)[1],
#                                                          quantile(all_data_clean$Area)[2])
#                                              & Area < mean(quantile(all_data_clean$Area)[4],
#                                                            quantile(all_data_clean$Area)[5])) &
#                               (Height > mean(quantile(all_data_clean$Height)[1],
#                                              quantile(all_data_clean$Height)[2])
#                                & Height < mean(quantile(all_data_clean$Height)[4],
#                                                quantile(all_data_clean$Height)[5])))
# }

filter_quantile <- subset(all_data_clean, (Area > 530000 & Area < 1050000) &
                                          (Height > 10000 & Height < 78000))

quantile(filter_quantile$Area)
quantile(filter_quantile$Height)
hist(filter_quantile$Area)
hist(filter_quantile$Height)

# Number of unique compounds after filtering
length(unique(all_data_clean$Compound))
length(unique(filter_quantile$Compound)) 


# Shapiro-Wilk Normality Test ---------------------------------------------
# Peak Area

# Peak Height


# Data Normalization -----------------------------------------------------

# # Pivot wider for PCA (may be unimportant) -----------------------------------------------------
# df_list[[1]] %>% 
#   pivot_longer(-c(Compound_and_sample, sample_name, Compound)) %>%
#   pivot_wider(names_from = Compound, values_from = value)

# i <- 1
# df_list_wider <- list()
# df_wider$sample_name <- file_list[[i]]
# df_wider <- relocate(df_wider, sample_name)
# df_list_wider[[i]] <- df_wider
# i <- i + 1


# PCA --------------------------------------------------------------------

# Change Compound column to row names
filter_quantile <- filter_quantile %>%
  column_to_rownames(., var = "Compound")

filter_quantile_pca <- PCA(filter_quantile[c(3:6)], scale.unit = TRUE, ncp = 5, graph = FALSE)

fviz_eig(filter_quantile_pca,
         addlabels = TRUE) 

var <- get_pca_var(filter_quantile_pca)
corrplot(var$cos2, is.corr = FALSE)
fviz_cos2(filter_quantile_pca, choice = "var", axes = 1:2)

ind <- get_pca_ind(filter_quantile_pca)
fviz_cos2(filter_quantile_pca, choice = "ind")

# Top 20 compounds with highest contribution
fviz_contrib(filter_quantile_pca, choice = "ind", axes = 1:2, top = 20) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))

# Investigate the significant compounds

