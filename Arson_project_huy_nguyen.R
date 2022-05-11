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

gasoline_df_list <- map(gasoline_file_list, read_excel, sheet = 'Results') 
# df_list <- map(file_list, read_excel, sheet = 'Results')
# Combine all data
# all_data <- bind_rows(df_list)


# Checkpoint for dimethylbenzene
# str_which(all_data$Compound, "(?=.*dimethyl)(?=.*benzene)") #2,4-Dinitro-1,3-dimethyl-benzene or (1,4-Dimethylpent-2-enyl)benzene


# Filtering out column bleed and solvent --------------------------------------

gasoline_df_list_clean <- map(gasoline_df_list, filtering, filter_list = c("^Carbon disulfide$", #df_list_clean
                                                                 "^Benzene$", 
                                                                 "Cyclotrisiloxane..hexamethyl",
                                                                 "Cyclotetrasiloxane..octamethyl",
                                                                 "^Toluene$",
                                                                 "^Ethylbenzene$",
                                                                 "Xylene")) 

# all_data_clean <- bind_rows(df_list_clean)

# Iterative loop sub-setting data based on cumulative sum(Percent_Height): cumulative sum(Percent_Height) < cumulative sum(Percent_Area)

# subset_df_list <- list()
gasoline_subset_df_list <- list() # subset_df_list

for (i in 1:length(gasoline_df_list_clean)) { # df_list_clean
  subset <- gasoline_df_list_clean[[i]]
  
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
    # slice data based on condition of cumulative sum of percent_height, limit ~ 80%
    if (sum(subset_clean[1:row_num,]$Percent_Height) > 0.95) {
      new_subset_clean <- slice_head(subset_clean, n = row_num)
      break
    }
  }
  # Assign new slice df to subset_df_list 
  gasoline_subset_df_list[[i]] <- new_subset_clean
}


# Export gasoline subset to csv file ------------------------------------------------------------------------------

# for (i in 1:length(gasoline_subset_df_list)) {
#   write.csv(gasoline_subset_df_list[[i]], file = paste0("C:/Users/huyng/Desktop/", 
#                                                         gasoline_file_list[[i]], ".csv"))
# }


# Add sample_name column to each subset_df

for (i in 1:length(gasoline_subset_df_list)) {
  gasoline_subset_df_list[[i]] <- gasoline_subset_df_list[[i]] %>%
    group_by(Compound) %>% 
    mutate(sample_name = gasoline_file_list[[i]])
}

# Combine all subset_df together

all_gasoline_subset_clean <- bind_rows(gasoline_subset_df_list)

# Identifying compound unique to each sample type -----------------------------------------------------------------

# all_subset_clean <- all_subset_clean %>%
#   add_column(sample_type = NA)
# 
# # Create factor column for 4 sample type: Gasoline, Diesel, GasolineComposite, DieselComposite
# for (sample in unique(all_subset_clean$sample_name)) {
#   if (str_detect(sample, "GasComp")) {
#     all_subset_clean[str_which(all_subset_clean$sample_name, 
#                                sample),]$sample_type <- "GasComp"
#   }
#   if (str_detect(sample, "DieselComp")) {
#     all_subset_clean[str_which(all_subset_clean$sample_name, 
#                                sample),]$sample_type <- "DieselComp"
#   }
#   if (str_ends(sample, "D.xlsx")) {
#     all_subset_clean[str_which(all_subset_clean$sample_name, 
#                                sample),]$sample_type <- "Diesel"
#   }
#   else {
#     all_subset_clean[str_which(all_subset_clean$sample_name, 
#                                sample),]$sample_type <- "Gas"
#   }
# }

# Conditional loop for sorting compound unique to each sample type
# Gasoline
# gasoline_unique <- data.frame(matrix(ncol = 11, nrow = 0))
# colnames(gasoline_unique) <- colnames(all_subset_clean)
# gasoline_unique$Compound <- as.character(gasoline_unique$Compound)
# gasoline_unique$sample_name <- as.character(gasoline_unique$sample_name)
# gasoline_unique$sample_type <- as.character(gasoline_unique$sample_type) #necessary????
# 
# for (compound in unique(all_subset_clean$Compound)) {
#   slice_df <- all_subset_clean[str_which(all_subset_clean$Compound, fixed(compound)),]
# 
#   if (unique(slice_df$sample_type) == "Gas") {
# 
#     gasoline_unique <- bind_rows(gasoline_unique, slice_df, .id = NULL)
#   }
# }
# 
# # Diesel
# diesel_unique <- data.frame(matrix(ncol = 10, nrow = 0))
# colnames(diesel_unique) <- colnames(all_subset_clean)
# diesel_unique$Compound <- as.character(diesel_unique$Compound)
# diesel_unique$sample_name <- as.character(diesel_unique$sample_name)
# diesel_unique$sample_type <- as.character(diesel_unique$sample_type)
# 
# for (compound in unique(all_subset_clean$Compound)) {
#   slice_df <- all_subset_clean[str_which(all_subset_clean$Compound, fixed(compound)),]
#   if (unique(slice_df$sample_name) == "Diesel") {
#     gasoline_unique <- bind_rows(gasoline_unique, slice_df, .id = NULL)
#   }
# }
# 
# # GasComp
# diesel_unique <- data.frame(matrix(ncol = 10, nrow = 0))
# colnames(diesel_unique) <- colnames(all_subset_clean)
# diesel_unique$Compound <- as.character(diesel_unique$Compound)
# diesel_unique$sample_name <- as.character(diesel_unique$sample_name)
# diesel_unique$sample_type <- as.character(diesel_unique$sample_type)
# 
# for (compound in unique(all_subset_clean$Compound)) {
#   slice_df <- all_subset_clean[str_which(all_subset_clean$Compound, fixed(compound)),]
#   if (unique(slice_df$sample_name) == "GasComp") {
#     gasoline_unique <- bind_rows(gasoline_unique, slice_df, .id = NULL)
#   }
# }
# 
# # DieselComp
# diesel_unique <- data.frame(matrix(ncol = 10, nrow = 0))
# colnames(diesel_unique) <- colnames(all_subset_clean)
# diesel_unique$Compound <- as.character(diesel_unique$Compound)
# diesel_unique$sample_name <- as.character(diesel_unique$sample_name)
# diesel_unique$sample_type <- as.character(diesel_unique$sample_type)
# 
# for (compound in unique(all_subset_clean$Compound)) {
#   slice_df <- all_subset_clean[str_which(all_subset_clean$Compound, fixed(compound)),]
#   if (unique(slice_df$sample_name) == "DieselComp") {
#     gasoline_unique <- bind_rows(gasoline_unique, slice_df, .id = NULL)
#   }
# }



# Dominant Compounds (high % area & height) found across all 39 samples ------------------------------------------------

all_gasoline_subset_clean_1st_filter <- data.frame(matrix(ncol = ncol(all_gasoline_subset_clean), nrow = 0))
colnames(all_gasoline_subset_clean_1st_filter) <- colnames(all_gasoline_subset_clean)
all_gasoline_subset_clean_1st_filter$Compound <- as.character(all_gasoline_subset_clean_1st_filter$Compound)
all_gasoline_subset_clean_1st_filter$sample_name <- as.character(all_gasoline_subset_clean_1st_filter$sample_name)
# all_subset_clean_1st_filter$sample_type <- as.character(all_subset_clean_1st_filter$sample_type)

for (compound_name in unique(all_gasoline_subset_clean$Compound)) { #all_subset_clean
  # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64
  if (sum(str_ends(compound_name, fixed(c("Citronellic acid", #0220F00191
                                          "4a,8,8-trimethyloctahydrocyclopropa(d)naphthalen-2(3H)-one" #0220F00694_check
                                          )))) > 0) {
    print(paste0("weird compound found ", compound_name))
    next
  }
  if (sum(str_starts(compound_name, fixed(c("1,4-Methanoazulen-9-one, decahydro-1,5,5,8a-tetramethyl",#0220F00694_check
                                            "Cyclohexanemethanol, 4-ethenyl",#0220F00694_check
                                            "(+-)-2-cis-4-trans-Abscisic acid", #0220F00694_check
                                            "1H-3a,7-Methanoazulene, octahydro-3,8,8-trimethyl-6-methylene", #0220F00694_check
                                            "2H-Pyran, tetrahydro-2-[2-(methylenecyclopropyl)ethoxy", #0220F00791
                                            "4-Octene, 2,6-dimethyl", # 0220F00591
                                            "Bicyclo[4.1.0]heptan-3-one, 4,7,7-trimethyl", #0220F00589 and 0220F00694_check
                                            "Bicyclo[4.1.0]heptane, 3,7,7-trimethyl", #0220F00487_check and 0220F00587
                                            "Cyclohexane, [6-cyclopentyl-3", #0220F00587 and 0220F00694_check
                                            "4H-1,3-Benzodioxin-4-one, hexahydro-4a,5-dimethyl", # 0220F00587
                                            "2H-Pyran, tetrahydro-2-[(tetrahydro-2-furanyl)methoxy", # 0220F00791
                                            "Nickel", #0220F00889
                                            "Tropinone", #0220F00889
                                            "1-Heptanol, 2,4-dimethyl", #0220F00989
                                            "Naphthalene, 1,1'-(1,10-decanediyl)bis[decahydro", #0220F00694_check
                                            "Naphthalene, decahydro-1,4a-dimethyl-7-(1-methylethyl", #0220F00694_check
                                            "Cyclopentane, 1,1'-[3-(2-cyclopentylethyl)-1,5-pentanediyl", # 0220F00694_check
                                            "1,4-Methanoazulen-9-one, decahydro-1,5,5,8a-tetramethyl", #0220F00694_check
                                            "1H-3a,7-Methanoazulene, 2,3,4,7,8,8a-hexahydro-3,6,8,8-tetramethyl", #0220F00694_check
                                            "Naphthalene, 1,2,3,5,6,7,8,8a-octahydro-1,8a-dimethyl-7-(1-methylethenyl", #0220F00694_check
                                            "Spiro[bicyclo[3.1.1]heptane-2,2'-oxirane], 6,6-dimethyl", #0220F00694_check
                                            "3a,7-Methano-3aH-cyclopentacyclooctene, decahydro-1,1,7-trimethyl", #0220F00694_check
                                            "Spiro[tricyclo[4.4.0.0(5,9)]decane-10,2'-oxirane], 1-methyl-4-isopropyl-7,8-dihydroxy", #0220F00694_check
                                            "1H-Cycloprop[e]azulene, 1a,2,3,5,6,7,7a,7b-octahydro-1,1,4,7-tetramethyl", #0220F00694_check
                                            "1,2,4-Metheno-1H-indene, octahydro-1,7a-dimethyl-5-(1-methylethyl", #0220F00694_check
                                            "Tricyclo[6.3.1.0(2,5)-]dodecan-1-ol, 4,4,8-trimethyl-, acetate", #0220F00694_check
                                            "1H-Cyclopropa[a]naphthalene, 1a,2,6,7,7a,7b-hexahydro-1,1,7,7a-tetramethyl", #0220F00694_check
                                            "(+)-Cycloisolongifol-5-ol", #0220F00694_check
                                            "Cyclopentane, 1,1'-[4-(3-cyclopentylpropyl)-1,7-heptanediyl]bis", #0220F00694_check
                                            "2aS,3aR,5aS,9bR)-2a,5a,9-Trimethyl-2a,4,5,5a,6,7,8,9b-octahydro-2H-naphtho[1,2-b]oxireno[2,3-c]furan", #0220F00694_check
                                            "Benzene, 1,1'-ethylidenebis[3,4-dimethyl" #0220F00694_check
                                            )))) > 0) {
    print(paste0("weird compound found ", compound_name))
    next
  }
  # First filter extract row index at the match and examine sample_name at that row index
  slice_df <- all_gasoline_subset_clean[str_which(all_gasoline_subset_clean$Compound, paste0("^", compound_name, "$")),] # DEBUG THIS! paste0("^", compound, "$")
  # !! Problem with Error in stri_detect_regex:  Syntax error in regex pattern
  # !! Problem with Error in stri_detect_regex: In a character range [x-y], x is greater than y. [name is too long]
  
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df$sample_name)) > 21) {
    # append the remove slice_df to all_subset_clean_similar
    all_gasoline_subset_clean_1st_filter <- bind_rows(all_gasoline_subset_clean_1st_filter, slice_df, .id = NULL)
  }
}


# examine the cumulative peak height and peak area per sample of compounds found across 39 samples
all_gasoline_subset_clean_1st_filter %>%
  group_by(sample_name) %>%
  summarise(cumulative_area = sum(Percent_Area))




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

