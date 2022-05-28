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
  # .[str_detect(., "D.xlsx")] %>%
  # .[str_detect(., "DieselComp")] %>%
  # .[str_detect(., "GasComp")] %>%
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


# Label compound type based on chemical structure -----------------------------------------------------------------
# REGEX REFERENCES:
# https://regenerativetoday.com/a-complete-beginners-guide-to-regular-expressions-in-r/
# https://en.wikibooks.org/wiki/R_Programming/Text_Processing#Regular_Expressions
# https://towardsdatascience.com/a-gentle-introduction-to-regular-expressions-with-r-df5e897ca432

# test1 <- df_list_clean[[1]] %>%
#   select(-ends_with(c("Area %", "Ion 1", "Ion 2", "Ion 3"))) %>%
#   mutate(compound_type = ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo", # ignore.case -> case insensitive
#                                 ifelse(grepl("bromo", Compound),"bromo",
#                                        ifelse(grepl("chloro", Compound, ignore.case = TRUE),"chloro",
#                                               ifelse(grepl("phospho", Compound),"phospho",
#                                                      ifelse(grepl("sulf", Compound),"sulfur",
#                                                             ifelse(grepl("amin", Compound),"amine",
#                                                                    ifelse(grepl("naphthal", Compound),"naphthalene",
#                                                                           ifelse(grepl("Benze", Compound), "benzene", "others")))))))))

# Iterative loop sub-setting data based on cumulative sum(Percent_Height): cumulative sum(Percent_Height) < cumulative sum(Percent_Area)

# subset_df_list <- list()
slice_df_list <- list() # subset_df_list

system.time({for (i in 1:length(df_list_clean)) { 
  df <- df_list_clean[[i]] %>%
    mutate(Percent_Area = Area/sum(Area)) %>%
    mutate(Percent_Height = Height/sum(Height)) %>%
    arrange(desc(Percent_Height), desc(Percent_Area)) # %>%
    # Compound column convert to rownames
    # mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height))
  
  
  # # subset data based on the largest number of iteration
  for (row_num in 1:nrow(df)) {
    # slice data based on condition of cumulative sum of percent_height, limit ~ 80%
    if (sum(df[1:row_num,]$Percent_Height) > 0.99) {
      new_df <- slice_head(df, n = row_num)
      break
    }
  }
  
  # Assign new slice df to subset_df_list 
  slice_df_list[[i]] <- new_df
  
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


# Dominant Compounds (high % area & height) found across samples ------------------------------------------------
# Similar compound from IL samples
all_similar_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0)) #all_gasoline_subset_clean_1st_filter
colnames(all_similar_compounds) <- colnames(all_subset_clean)
all_similar_compounds$Compound <- as.character(all_similar_compounds$Compound)
all_similar_compounds$sample_name <- as.character(all_similar_compounds$sample_name)
# all_subset_clean_1st_filter$sample_type <- as.character(all_subset_clean_1st_filter$sample_type)

# Distinctive compound for IL samples
all_different_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0)) #unique_gasoline_subset
colnames(all_different_compounds) <- colnames(all_subset_clean)
all_different_compounds$Compound <- as.character(all_different_compounds$Compound)
all_different_compounds$sample_name <- as.character(all_different_compounds$sample_name)

system.time({for (compound_name in unique(all_subset_clean$Compound)) { #all_subset_clean
  # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64

  # Filter 1: baseR::grepl 
  slice_df1 <- all_subset_clean[which(grepl(compound_name, all_subset_clean$Compound, fixed = TRUE)),]
  
  # !!! --> If try() produces error, then next iteration in for loop
  if (class(try(which(grepl(paste0("^", compound_name, "$"), all_subset_clean$Compound)), silent = TRUE)) %in% "try-error") {
    # put weird compound names in unique subset
    all_different_compounds <- bind_rows(all_different_compounds, slice_df1, .id = NULL)
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
    all_similar_compounds <- bind_rows(all_similar_compounds, slice_df2, .id = NULL)
  }
  else {
    all_different_compounds <- bind_rows(all_different_compounds, slice_df2, .id = NULL)
  }
  rm(slice_df2)
}})

length(unique(all_similar_compounds$Compound))
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

length(unique(all_different_compounds$Compound))
# 963 compounds unique for 5 diesel compounds

# examine the cumulative peak height and peak area per sample of compounds found across 39 samples


# Filtering unique compound for each sample -----------------------------------------------------------------------
all_unique_compounds <- data.frame(matrix(ncol = ncol(all_different_compounds), nrow = 0))
colnames(all_unique_compounds) <- colnames(all_different_compounds)
all_unique_compounds$Compound <- as.character(all_unique_compounds$Compound)
all_unique_compounds$sample_name <- as.character(all_unique_compounds$sample_name)

system.time({
  for (compound_name in unique(all_different_compounds$Compound)) {
    slice_df <- all_different_compounds[which(grepl(compound_name, all_different_compounds$Compound, fixed = TRUE)),]
    if (length(unique(slice_df$sample_name)) < 2) {
      # append the remove slice_df to all_subset_clean_similar
      all_unique_compounds <- bind_rows(all_unique_compounds, slice_df, .id = NULL)
    } 
  }
}) 

length(unique(all_unique_compounds$Compound))
# 832 compounds unique for  gasoline compounds

# PCA on indi_IL_type ---------------------------------------------------------------------------------------------

unique_subsetdf <- all_unique_compounds %>%                                    # ignore.case -> case insensitive
  # mutate(compound_type = ifelse(grepl("bromo", Compound, ignore.case = TRUE),"bromo",
  #                               ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo",
  #                                      ifelse(grepl("chlor", Compound, ignore.case = TRUE),"chloro",
  #                                             ifelse(grepl("phosph", Compound, ignore.case = TRUE),"phospho",
  #                                                    ifelse(grepl("sulf", Compound, ignore.case = TRUE),"sulfur",
  #                                                           ifelse(grepl("amin", Compound, ignore.case = TRUE),"amine",
  #                                                                  ifelse(grepl("naphthal", Compound, ignore.case = TRUE),"naphthalene",
  #                                                                         ifelse(grepl("Benze", Compound, ignore.case = TRUE), "benzene","others"))))))))) %>%
  # filter(!grepl("others", compound_type)) %>%
  mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height)) %>%
  column_to_rownames(., var = "Compound")

# MUST CONVERT GROUPING COLUMN to factor type, otherwise error "undefined columns selected" will happen for 'habillage'
unique_subsetdf$sample_name <- factor(unique_subsetdf$sample_name, levels = c(unique(unique_subsetdf$sample_name)))

# Histogram plot for Percent Height 

hist(unique_subsetdf$Percent_Height, xlim = c(0, 0.002), breaks = 3000)
quantile(unique_subsetdf$Percent_Height)
filter_quantile <- subset(unique_subsetdf, 
                          Percent_Height >  0.00021 &
                            Percent_Height < 0.002)
length(unique(filter_quantile$sample_name))

# Manual checkpoint for correlation of different variables (RT1,Rt2, %Area, %height, etc.) via covariance matrix 
# CAUTION!!! - If the variables are not strongly correlated (abs. covariance value must > 0.75), then there is no point to use PCA
view(cov(scale(filter_quantile[c(1:12)], center = TRUE, scale = TRUE)))
view(cor(scale(filter_quantile[c(1:12)], center = TRUE, scale = TRUE), method = "spearman"))

# Compare the correlation matrix our data to iris data
view(cor(scale(iris[c(1:4)], center = TRUE, scale = TRUE), method = "spearman"))

filter_quantile_pca <- PCA(filter_quantile[c(3,4,11,12)]
                           , scale.unit = TRUE, graph = FALSE)
# Investigate grouping of compounds- REFERENCE: https://pca4ds.github.io/biplot-and-pca.html
fviz_pca_biplot(filter_quantile_pca, repel = TRUE, label = "var",
                       habillage = filter_quantile$sample_name,
                       # addEllipses=TRUE,
                       dpi = 480)

# Compare the PCA of our data to iris data
iris_pca <- PCA(iris[c(1:4)], scale.unit = TRUE)
fviz_pca_biplot(iris_pca, repel = TRUE, label = "var",
                habillage = iris$Species,
                dpi = 480)
# ggsave(paste0(getwd(), "/PCA graphs/dieselcomp_pca.png"),
#        pca,
#        height = 8,
#        width = 15)

# PCA using only Percent_height value on compounds that exist in only one sample
subset_filterquantile <- all_similar_compounds %>%
  # mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", `Ion 1`, "-", `Ion 2`, "-", `Ion 3`)) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  group_by(sample_name, Compound) %>% 
  summarise(across(c(Percent_Area, Percent_Height), mean)) %>%
  # filter(Percent_Height >  0.00021 & Percent_Height < 0.002) %>%
  # select(Compound, Percent_Height, sample_name) %>%
  pivot_wider(names_from = Compound, values_from = c(Percent_Area, Percent_Height))

subset_filterquantilePCA <- PCA(subset_filterquantile[c(2:69)], scale.unit = TRUE, graph = FALSE)

fviz_eig(subset_filterquantilePCA,
         addlabels = TRUE)

fviz_pca_biplot(subset_filterquantilePCA, 
                # repel = TRUE,
                label = "ind",
                axes = c(1,2),
                habillage = subset_filterquantile$sample_name,
                # addEllipses=TRUE,
                dpi = 480)

# PCA --------------------------------------------------------------------

# Individual IL files
# Data frame for sorting percent_area & percent_height from highest to lowest
for (i in 1:length(df_list_clean)) {
  testdf <- df_list_clean[[i]] %>%
    mutate(Percent_Area = Area/sum(Area)) %>%
    mutate(Percent_Height = Height/sum(Height)) %>%
    arrange(desc(Percent_Height), desc(Percent_Area)) %>%
    # Compound column convert to rownames
    mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height))
  
  
  # # subset data based on the largest number of iteration
  for (row_num in 1:nrow(testdf)) {
    # slice data based on condition of cumulative sum of percent_height, limit ~ 80%
    if (sum(testdf[1:row_num,]$Percent_Height) > 0.99) {
      testdf <- slice_head(testdf, n = row_num)
      break
    }
  }
  
  # cov_df_scaled <- cov(scale(testdf %>%
  #                      column_to_rownames(., var = "Compound"),
  #                    center = TRUE, scale = TRUE))
  # view(cov_df_scaled)
  # Grouping compound types
  testdf <- all_subset_clean %>%                                    # ignore.case -> case insensitive
    # mutate(compound_type = ifelse(grepl("bromo", Compound, ignore.case = TRUE),"bromo",
    #                               ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo",
    #                                      ifelse(grepl("chlor", Compound, ignore.case = TRUE),"chloro",
    #                                             ifelse(grepl("phosph", Compound, ignore.case = TRUE),"phospho",
    #                                                    ifelse(grepl("sulf", Compound, ignore.case = TRUE),"sulfur",
    #                                                           ifelse(grepl("amin", Compound, ignore.case = TRUE),"amine",
    #                                                                  ifelse(grepl("naphthal", Compound, ignore.case = TRUE),"naphthalene",
    #                                                                         ifelse(grepl("Benze", Compound, ignore.case = TRUE), "benzene","others"))))))))) %>%
    # filter(!grepl("others", compound_type)) %>%
    column_to_rownames(., var = "Compound")
  
  # MUST CONVERT GROUPING COLUMN to factor type, otherwise error "undefined columns selected" will happen for 'habillage'
  testdf$sample_name <- factor(testdf$sample_name, levels = c(unique(testdf$sample_name)))
  
  filter_quantile_pca <- PCA(testdf[c(3,4,11,12)], scale.unit = TRUE, graph = FALSE)
  
  # Investigate grouping of compounds
  # REFERENCE: https://pca4ds.github.io/biplot-and-pca.html
  ILR_pca <- fviz_pca_biplot(filter_quantile_pca, repel = TRUE, label = "var",
                  habillage = testdf$sample_name,
                  # palette = "Dark2",
                  addEllipses=TRUE,
                  dpi = 480)
  # Scree plot
  fviz_eig(filter_quantile_pca, addlabels = TRUE)
  
  # Cos2 aka. Quality of representation
  # Source: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials
  var <- get_pca_var(filter_quantile_pca)
  corrplot(var$cos2, is.corr = FALSE)
  
  fviz_cos2(filter_quantile_pca, choice = "var", axes = 1, top = 20)
  
  # Top variables (RT1, RT2,etc.) and compounds with highest contribution
  fviz_contrib(filter_quantile_pca, choice = "var", axes = 1) # contrib of var to PC1
  fviz_contrib(filter_quantile_pca, choice = "var", axes = 2) # contrib of var to PC2
  fviz_contrib(filter_quantile_pca, choice = "ind", axes = 2, top = 20) + # contrib of indi to PC2
    theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))
  # ggsave(paste0(getwd(), "/PCA graphs/ILR_pca.png"),
  #        ILR_pca,
  #        height = 8,
  #        width = 15)
}
