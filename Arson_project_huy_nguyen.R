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

# Pipe operator for list
gasoline_file_list <- file_list %>% 
  .[!str_ends(., "D.xlsx")] %>%
  .[!str_detect(., "DieselComp")] %>%
  .[!str_detect(., "GasComp")] #%>%
  # .[!str_ends(., "check.xlsx")]

# gasoline_df_list <- map(gasoline_file_list, read_excel, sheet = 'Results') 
df_list <- map(file_list, read_excel, sheet = 'Results')
# Combine all data
# all_data <- bind_rows(df_list)


# Checkpoint for dimethylbenzene
# str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene


# Filtering out column bleed and solvent --------------------------------------

df_list_clean <- map(df_list, filtering, filter_list = c("^Carbon disulfide$", #gasoline_df_list_clean, gasoline_df_list
                                                                 "^Benzene$", 
                                                                 "Cyclotrisiloxane..hexamethyl",
                                                                 "Cyclotetrasiloxane..octamethyl",
                                                                 "^Toluene$",
                                                                 "^Ethylbenzene$",
                                                                 "Xylene")) 

# all_data_clean <- bind_rows(df_list_clean)

# Iterative loop sub-setting data based on cumulative sum(Percent_Height): cumulative sum(Percent_Height) < cumulative sum(Percent_Area)

# subset_df_list <- list()
subset_df_list <- list() # subset_df_list

for (i in 1:length(df_list_clean)) { # df_list_clean
  subset <- df_list_clean[[i]]
  
  # Optional !! Remove unnecessary columns in data
  subset_clean <- subset %>%
    select(-ends_with(c("Area %", "Ion 1", "Ion 2", "Ion 3")))
  
  # METHOD 1: BASED ON CUMULATIVE SUM OF PERCENTAGE HEIGHT
  
  # Calculate percentage Peak Area and Peak Height
  subset_clean$Percent_Area <- subset_clean$Area/sum(subset_clean$Area)
  subset_clean$Percent_Height <- subset_clean$Height/sum(subset_clean$Height)

  # Data frame for sorting percent_area & percent_height from highest to lowest
  subset_clean <- subset_clean %>%
    arrange(desc(Percent_Height), desc(Percent_Area))

  # subset data based on the largest number of iteration
  for (row_num in 1:nrow(subset_clean)) {
    # slice data based on condition of cumulative sum of percent_height, limit ~ 80%
    if (sum(subset_clean[1:row_num,]$Percent_Height) > 0.95) {
      new_subset_clean <- slice_head(subset_clean, n = row_num)
      break
    }
  }
  # Assign new slice df to subset_df_list 
  subset_df_list[[i]] <- new_subset_clean
  
  # METHOD 2: BASED ON CUMULATIVE DISTRIBUTION OF AREA UNDER THE CURVE
}


# Export gasoline subset to csv file ------------------------------------------------------------------------------

# for (i in 1:length(gasoline_subset_df_list)) {
#   write.csv(gasoline_subset_df_list[[i]], file = paste0("C:/Users/huyng/Desktop/", 
#                                                         gasoline_file_list[[i]], ".csv"))
# }


# Add sample_name column to each subset_df

for (i in 1:length(subset_df_list)) {
  subset_df_list[[i]] <- subset_df_list[[i]] %>%
    group_by(Compound) %>% 
    mutate(sample_name = file_list[[i]]) #gasoline_file_list
}

# Combine all subset_df together

all_subset_clean <- bind_rows(subset_df_list) # all_subset_clean, all_gasoline_subset_clean, gasoline_subset_df_list

# Dominant Compounds (high % area & height) found across all 39 samples ------------------------------------------------
# Similar compound from IL samples
all_subset_clean_1st_filter <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0)) #all_gasoline_subset_clean_1st_filter
colnames(all_subset_clean_1st_filter) <- colnames(all_subset_clean)
all_subset_clean_1st_filter$Compound <- as.character(all_subset_clean_1st_filter$Compound)
all_subset_clean_1st_filter$sample_name <- as.character(all_subset_clean_1st_filter$sample_name)
# all_subset_clean_1st_filter$sample_type <- as.character(all_subset_clean_1st_filter$sample_type)

# Distinctive compound for IL samples
unique_subset <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0)) #unique_gasoline_subset
colnames(unique_subset) <- colnames(all_subset_clean)
unique_subset$Compound <- as.character(unique_subset$Compound)
unique_subset$sample_name <- as.character(unique_subset$sample_name)

for (compound_name in unique(all_subset_clean$Compound)) { #all_subset_clean
  # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64

  # Filter 1: baseR::grepl 
  slice_df1 <- all_subset_clean[which(grepl(compound_name, all_subset_clean$Compound, fixed = TRUE)),]
  
  # !!! --> If try() produces error, then next iteration in for loop
  if (class(try(str_which(slice_df1$Compound, paste0("^", compound_name, "$")), silent = TRUE)) %in% "try-error") {
    unique_subset <- bind_rows(unique_subset, slice_df1, .id = NULL)
    next
  }
  # Filter 2: stringr::str_which
  else {
    slice_df2 <- slice_df1[str_which(slice_df1$Compound, paste0("^", compound_name, "$")),]
  }
  
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df2$sample_name)) > 15) {
    # append the remove slice_df to all_subset_clean_similar
    all_subset_clean_1st_filter <- bind_rows(all_subset_clean_1st_filter, slice_df2, .id = NULL)
  }
  else {
    unique_subset <- bind_rows(unique_subset, slice_df2, .id = NULL)
  }
}

length(unique(all_subset_clean_1st_filter$Compound)) 
# more than 15 out of 28 gasoline samples share 96 common compounds
# more than 15 out of 39 IL samples share 150 common compounds
unique_subset <- unique_subset %>%
  arrange(desc(sample_name), desc(MF)) #%>%
  #arrange(desc(Percent_Area), desc(Percent_Height))
length(unique(unique_subset$sample_name))
# examine the cumulative peak height and peak area per sample of compounds found across 39 samples
# all_gasoline_subset_clean_1st_filter %>%
#   group_by(sample_name) %>%
#   summarise(cumulative_area = sum(Percent_Area))

unique_subset_comp_freq <- unique_subset %>% group_by(Compound) %>% summarise(freq_occur = frequency(Compound))
max(unique_subset_comp_freq$freq_occur, na.rm = TRUE)
# Filtering unique compound for each sample -----------------------------------------------------------------------




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

