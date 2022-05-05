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
library(stringi)

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
 # Notin function
`%notin%` <- Negate(`%in%`)

# Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')
df_list <- map(file_list, read_excel, sheet = 'Results')

# Combine all data
all_data <- bind_rows(df_list)

# Checkpoint for dimethylbenzene
# str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene


# Filtering out column bleed and solvent --------------------------------------

df_list_clean <- map(df_list, filtering, filter_list = c("Carbon.disulfide",
                                        "Benzene - ", 
                                        "Cyclotrisiloxane..hexamethyl",
                                        "Cyclotetrasiloxane..octamethyl",
                                        "Toluene",
                                        "Ethylbenzene",
                                        "Xylene")) 

all_data_clean <- bind_rows(df_list_clean)

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
  # Assign new slice df to subset_df_list 
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

# filtering similar compound across all 39 samples in subset_df_list
all_subset_clean_similar <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(all_subset_clean_similar) <- colnames(all_subset_clean)
all_subset_clean_similar$Compound <- as.character(all_subset_clean_similar$Compound)
all_subset_clean_similar$sample_name <- as.character(all_subset_clean_similar$sample_name)

for (compound in unique(all_subset_clean$Compound)) {
  # ignore long compound names ????
  if (nchar(compound) > 48) {
    next
  }
  
  compound <- paste0("^",compound,"$")
  # extract row index at the match and examine sample_name at that row index
  slice_df <- all_subset_clean[str_which(all_subset_clean$Compound, compound,] # DEBUG THIS!
  # Problem with Error in stri_detect_regex:  Syntax error in regex pattern
  # Problem with Error in stri_detect_regex: In a character range [x-y], x is greater than y. [name is too long]
  
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df$sample_name)) == 39) {
    # append the remove slice_df to all_subset_clean_similar
    all_subset_clean_similar <- bind_rows(all_subset_clean_similar, slice_df, .id = NULL)
    print(dim(all_subset_clean_similar))
  }
}
 


# Check number of unique compound 
length(unique(all_subset_clean_similar$Compound))

# examine the cumulative peak height and peak area per sample of compounds found across 39 samples
all_subset_clean_similar %>%
  group_by(sample_name) %>%
  summarise(cumulative_area = sum(Percent_Area)) 
# Error!!! cumulative_area value of ecah sample here should be less than sum percent_area 
# of each sample in subset_df_list



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

