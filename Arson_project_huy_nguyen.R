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
library(purrr)
library(stringr)

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

# all_data_clean <- bind_rows(df_list_clean)

# Iterative loop subsetting data based on cumulative sum of Percent_Height

subset_df_list <- list()

for (i in 1:length(df_list_clean)) {
  subset <- df_list_clean[[i]]
  
  # Optional !! Remove unnecessary columns in data
  subset_clean <- subset %>%
    select(-ends_with(c("Area %", "Ion 1", "Ion 2", "Ion 3")))

  # Calculate percentage Peak Area and Peak Height
  subset_clean$Percent_Area <- subset_clean$Area/sum(subset_clean$Area)
  subset_clean$Percent_Height <- subset_clean$Height/sum(subset_clean$Height)

  # Data frame for sorting percent_area & percent_height from highest to lowest
  subset_clean <- subset_clean %>%
    arrange(desc(Percent_Height), desc(Percent_Area))

  # subset data based on the largest number of iteration
  for (row_num in 1:nrow(subset_clean)) {
    # slice data based on condition of cumulative sum of percent_height, limit = 80%
    if (sum(subset_clean[1:row_num,]$Percent_Height) > 0.8) {
      new_subset_clean <- slice_head(subset_clean, n = row_num)
      break
    }
  }
  # Assign the 
  subset_df_list[[i]] <- new_subset_clean
}

# Add sample_name column to each subset_df

for (i in 1:length(subset_df_list)) {
  subset_df_list[[i]] <- subset_df_list[[i]] %>%
    group_by(Compound) %>% 
    mutate(sample_name = file_list[[i]])
}

# Combine all subset_df together

all_subset_clean <- bind_rows(subset_df_list)

# filtering similar compound across all 39 files in subset_df_list
all_subset_clean_similar <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(all_subset_clean_similar) <- colnames(all_subset_clean)
all_subset_clean_similar$Compound <- as.character(all_subset_clean_similar$Compound)
all_subset_clean_similar$sample_name <- as.character(all_subset_clean_similar$sample_name)

for (compound in unique(all_subset_clean$Compound)) {

  # extract row index at the match and examine sample_name at that row index
  slice_df <- all_subset_clean[str_detect(all_subset_clean$Compound, coll(compound)),]
  
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df$sample_name)) == 39) {
    all_subset_clean_similar <- bind_rows(all_subset_clean_similar, slice_df, .id = NULL)
  }
}

# Check number of unique compound 
length(unique(all_subset_clean_similar$Compound))

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

