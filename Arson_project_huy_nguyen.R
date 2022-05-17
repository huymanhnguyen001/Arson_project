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

options(ggrepel.max.overlaps = 20)
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

# Pipe operator for isolating IL types
indi_IL_file_list <- file_list %>%
  # .[!str_detect(., "D.xlsx")] %>%
  # .[!str_detect(., "DieselComp")] %>%
  # .[!str_detect(., "GasComp")] %>%
  .[!str_ends(., "check.xlsx")]

# Import IL samples to list
df_list <- map(indi_IL_file_list, read_excel, sheet = 'Results')

# Combine all data
# all_data <- bind_rows(df_list)


# Checkpoint for dimethylbenzene
# str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene


# Filtering out column bleed and solvent --------------------------------------

df_list_clean <- map(df_list, filtering, filter_list = c("^Carbon disulfide$", 
                                                                 "^Benzene$", 
                                                                 "Cyclotrisiloxane..hexamethyl",
                                                                 "Cyclotetrasiloxane..octamethyl",
                                                                 "^Toluene$",
                                                                 "^Ethylbenzene$",
                                                                 "Xylene")) 

# all_data_clean <- bind_rows(df_list_clean)

# Iterative loop sub-setting data based on cumulative sum(Percent_Height): cumulative sum(Percent_Height) < cumulative sum(Percent_Area)

# subset_df_list <- list()
slice_df_list <- list() # subset_df_list

system.time({for (i in 1:length(df_list_clean)) { 
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
    if (sum(subset_clean[1:row_num,]$Percent_Height) > 0.99) {
      new_subset_clean <- slice_head(subset_clean, n = row_num)
      break
    }
  }
  rm(subset_clean)
  # Assign new slice df to subset_df_list 
  slice_df_list[[i]] <- new_subset_clean
  
  rm(new_subset_clean)
  
  # METHOD 2: BASED ON CUMULATIVE DISTRIBUTION OF AREA UNDER THE CURVE
}})

# Add sample_name column to each subset_df

for (i in 1:length(slice_df_list)) {
  slice_df_list[[i]] <- slice_df_list[[i]] %>%
    group_by(Compound) %>% 
    mutate(sample_name = indi_IL_file_list[[i]])
}

# Export gasoline subset to csv file ------------------------------------------------------------------------------

# for (i in 1:length(subset_df_list)) {
#   write.csv(subset_df_list[[i]], file = paste0("C:/Users/huyng/Desktop/",
#                                                         diesel_file_list[[i]], ".csv"))
# }

# Combine all subset_df together

all_subset_clean <- bind_rows(slice_df_list) 


# Create unique observation name with Compound + Sample_name + MF + RMF -----------------------------------------------
for (i in 1:length(df_list)) {
  df_list[[i]] <- df_list[[i]] %>% 
    group_by(Compound) %>% 
    summarise(across(everything(), mean))#  %>%
  # mutate(sample_name = file_list[[i]]) %>%
  # relocate(sample_name) %>%
  # mutate(Compound_and_sample = paste(Compound, "-", sample_name)) %>%
  # relocate(Compound_and_sample)
}


# Dominant Compounds (high % area & height) found across samples ------------------------------------------------
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

system.time({for (compound_name in unique(all_subset_clean$Compound)) { #all_subset_clean
  # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64

  # Filter 1: baseR::grepl 
  slice_df1 <- all_subset_clean[which(grepl(compound_name, all_subset_clean$Compound, fixed = TRUE)),]
  
  # !!! --> If try() produces error, then next iteration in for loop
  if (class(try(which(grepl(paste0("^", compound_name, "$"), all_subset_clean$Compound)), silent = TRUE)) %in% "try-error") {
    # put weird compound names in unique subset
    unique_subset <- bind_rows(unique_subset, slice_df1, .id = NULL)
    next
  } 
  # Filter 2: stringr::str_which
  else {
    slice_df2 <- slice_df1[which(grepl(paste0("^", compound_name, "$"), slice_df1$Compound)),]
  }
  rm(slice_df1)
  rm(compound_name)
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df2$sample_name)) > (length(indi_IL_file_list) - 1)) {
    # append the remove slice_df to all_subset_clean_similar
    all_subset_clean_1st_filter <- bind_rows(all_subset_clean_1st_filter, slice_df2, .id = NULL)
  }
  else {
    unique_subset <- bind_rows(unique_subset, slice_df2, .id = NULL)
  }
  rm(slice_df2)
}})

length(unique(all_subset_clean_1st_filter$Compound))
# When include 99% of cumulative peak height, all diesel samples share 304 compounds in common
# When include 99% of cumulative peak height, all gasoline samples share 39 compounds in common
# When include 99% of cumulative peak height, all diesel composite samples share 357 compounds in common
# When include 99% of cumulative peak height, all gasoline composite samples share 248 compounds in common
# When include 99% of cumulative peak height, all IL samples share 29 compounds in common

#--------------------------------------------------------------------------------------
# 28 out of 28 gasoline samples share 25 common compounds
# more than 15 out of 28 gasoline samples share 96 common compounds
# 5 out of 5 diesel samples share 50 common compounds
# more than 15 out of 39 IL samples share 150 common compounds

length(unique(unique_subset$Compound))
# 963 compounds unique for 5 diesel compounds

# examine the cumulative peak height and peak area per sample of compounds found across 39 samples
# all_gasoline_subset_clean_1st_filter %>%
#   group_by(sample_name) %>%
#   summarise(cumulative_area = sum(Percent_Area))
# 
# unique_subset_comp_freq <- unique_subset %>% group_by(Compound) %>% summarise(freq_occur = frequency(Compound))
# max(unique_subset_comp_freq$freq_occur, na.rm = TRUE)

# Filtering unique compound for each sample -----------------------------------------------------------------------
unique_subset_clean <- data.frame(matrix(ncol = ncol(unique_subset), nrow = 0))
colnames(unique_subset_clean) <- colnames(unique_subset)
unique_subset_clean$Compound <- as.character(unique_subset_clean$Compound)
unique_subset_clean$sample_name <- as.character(unique_subset_clean$sample_name)

system.time({
  for (compound_name in unique(unique_subset$Compound)) {
    slice_df <- unique_subset[which(grepl(compound_name, unique_subset$Compound, fixed = TRUE)),]
    if (length(unique(slice_df$sample_name)) < 2) {
      # append the remove slice_df to all_subset_clean_similar
      unique_subset_clean <- bind_rows(unique_subset_clean, slice_df, .id = NULL)
    } 
  }
}) 

length(unique(unique_subset_clean$Compound))
# 832 compounds unique for  gasoline compounds


# PCA --------------------------------------------------------------------

# Change Compound column to row names - maybe redundant
# filter_quantile <- all_subset_clean %>%
#   column_to_rownames(., var = "Compound")
slice_df31 <- slice_df_list[[31]] %>%
  mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height)) %>%
  column_to_rownames(., var = "Compound")

filter_quantile_pca <- PCA(slice_df31[c(3,4,7,8)], scale.unit = TRUE, ncp = 3, graph = TRUE)

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
slice_df31_biplot <- fviz_pca_biplot(filter_quantile_pca, 
                                    repel = TRUE, 
                                    col.var = "#2E9FDF",
                                    col.ind = "red")
                                    # habillage = slice_df1$Compound_type)
ggsave(paste0(getwd(), "/PCA graphs/GasComp samples/slice_df31_biplot.png"), 
       slice_df31_biplot,
       dpi = 240)
