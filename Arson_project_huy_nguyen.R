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
library(ChemoSpec) # ANOVA-PCA combo
library(DiscriMiner) # PLS-DA combo
library(missMDA)
library(tsne)
library(Rtsne)
library(rgl)
library(collapse)
library(plotly)
library(umap)
library(sqldf)
library(svMisc)
library(multiway)
library(foreach)
library(doSNOW)
library(writexl)
library(rapportools)
vignette("parallel")
options(ggrepel.max.overlaps = 300)
set.seed(12345)
cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)
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
filtering <- function(data, filter_list) { #, percentile, column_list
  clean_data <- copy(data)
  for (ele in filter_list) {
    clean_data <- clean_data %>%
      filter(!grepl(ele, Compound))
  }
  return(clean_data)
}

# Notin function
`%notin%` <- Negate(`%in%`)

# The InDaPCA function performs PCA with incomplete data -- WORSE THAN imputePCA
InDaPCA <- function(RAWDATA){
  
  #scaling
  
  X <- scale(RAWDATA, center = T, scale = T)
  
  #correlation
  
  C<-cor(X, use="pairwise.complete.obs")
  
  #Eigenvalue
  
  Eigenvalues<-eigen(C)$values
  
  Eigenvalues.pos<-Eigenvalues[Eigenvalues>0]
  
  Eigenvalues.pos.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
  
  #Eigenvectors
  
  V <- eigen(C)$vectors
  
  #Principal components
  
  X2<-X
  
  X2[is.na(X2)] <- 0
  
  PC <- as.matrix(X2) %*% V
  
  #object.standardized
  
  PCstand1 <- PC[,Eigenvalues>0]/sqrt(Eigenvalues.pos)[col(PC[,Eigenvalues>0])]
  
  PCstand2 <- PCstand1 / sqrt(nrow(PC) - 1)
  
  #loadings
  
  #loadings<-V%*%diag(sqrt(Eigenvalues.pos))
  
  loadings<-cor(X,PC,use="pairwise.complete.obs")
  
  #arrows for biplot
  
  arrows<-cor(X,PC,use="pairwise.complete.obs")*sqrt(nrow(X) - 1)
  
  #output
  
  PCA <- list()
  
  PCA$Correlation.matrix<-C
  
  PCA$Eigenvalues<-Eigenvalues
  
  PCA$Positive.Eigenvalues<-Eigenvalues.pos
  
  PCA$Positive.Eigenvalues.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
  
  PCA$Square.root.of.eigenvalues <- sqrt(Eigenvalues.pos)
  
  PCA$Eigenvectors<-V
  
  PCA$Component.scores<-PC
  
  PCA$Variable.scores<-loadings
  
  PCA$Biplot.objects<-PCstand2
  
  PCA$Biplot.variables<-arrows
  
  return(PCA)
  
}

# Grouping compounds based on RT1, RT2, and Ion1
grouping_comp <- function(data) {

  # create empty list, each sub-list is a compound group with following criteria:
  # RT1 threshold +- 0.4
  # RT2 threshold +- 0.1
  # Ion1 threshold +-0.1
  
  # Initialize the compound_group column filled with NA values
  data$compound_group <- NA
  i <- 1
  for (row in 1:nrow(data)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    rt1 <- data[row,]$RT1
    rt2 <- data[row,]$RT2
    ion1 <- data[row,]$`Ion 1`
    idx <- which(data$RT1 < (rt1 + 0.4) & data$RT1 > (rt1 - 0.4) & 
                 data$RT2 < (rt2 + 0.1) & data$RT2 > (rt2 - 0.1) & 
                 data$`Ion 1` < (ion1 + 0.1) & data$`Ion 1` > (ion1 - 0.1) & 
                 is.na(data$compound_group))
    if (identical(idx, integer(0))) {
      next
    }
    else {
      data[idx, "compound_group"] <- paste0("Compound_", i)
    }
    rm(rt1)
    rm(rt2)
    rm(ion1)
    i <- i + 1
  }
  return(data)
}

# Filtering similar and unique compound
comp_filter <- function(data, file_list) {
  all_similar_compounds_idx <- c()
  all_different_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$compound_group)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    
    idx <- which(grepl(paste0("^", comp_grp, "$"), data$compound_group))
    
    if (length(unique(data[idx,]$sample_name)) > (length(file_list) - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$sample_name)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_different_compounds_idx <- c(all_different_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_different_compounds_idx, all_unique_compounds_idx))
}

# Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')

# Pipe operator for isolating IL types
indi_IL_file_list <- file_list %>%
  # .[str_detect(., "D")]  %>%
  # .[str_detect(., "DieselComp")] %>%
  # .[str_detect(., "GasComp")] %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")]

# Import IL samples to list
df_list <- map(indi_IL_file_list, read_excel, sheet = "Results")

# Filtering out column bleed and solvent --------------------------------------
df_list_clean <- map(df_list, filtering, filter_list = c("^Carbon disulfide$", 
                                                                 "^Benzene$", 
                                                                 "Cyclotrisiloxane..hexamethyl",
                                                                 "Cyclotetrasiloxane..octamethyl",
                                                                 "^Toluene$",
                                                                 "^Ethylbenzene$",
                                                                 "Xylene")) 


# Iterative loop sub-setting data based on cumulative sum(Percent_Height)  --------
slice_df_list <- list() 

system.time({for (i in 1:length(df_list_clean)) { 
  df <- df_list_clean[[i]] %>%
    mutate(Percent_Area = Area/sum(Area)) %>%
    mutate(Percent_Height = Height/sum(Height)) %>%
    transmute(rowwise(.), Percent_Area = sort(c_across(Percent_Area), decreasing = TRUE)) %>%
    transmute(rowwise(.), Percent_Height = sort(c_across(Percent_Height), decreasing = TRUE))
  
  # subset data based on the largest number of iteration
  for (row_num in 1:nrow(df)) {
    # slice data based on condition of cumulative sum of percent_height
    if (sum(df[1:row_num,]$Percent_Height) > 0.99) {
      new_df <- slice_head(df, n = row_num)
      break
    }
  }
 
  # Add sample_name column 
  slice_df_list[[i]] <- new_df %>%
    mutate(sample_name = indi_IL_file_list[[i]]) %>%
    mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp",
                              ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                     ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
}})

# Grouping compounds based on RT1, RT2, Ion1 -----------------------------------------------------------------------
# Combine all subset_df together
all_subset_clean <- bind_rows(slice_df_list)

all_subset_clean_grouped <- grouping_comp(all_subset_clean)

# Export df for later use
write_xlsx(all_subset_clean_grouped, paste0(getwd(), "/grouping_compounds_200622.xlsx"))
# testing_import <- read_excel(paste0(getwd(), "/grouping_compounds.xlsx"))


# Similar Compounds (high % area & height) found across samples ------------------------------------------------
# Approach 1: Using compound "groups" by RT1, RT2, Ion1 Threshold
all_similar_compounds_idx1 <- c()
all_different_compounds_idx1 <- c()
all_unique_compounds_idx1 <- c()

system.time({for (comp_grp in unique(all_subset_clean_grouped$compound_group)) {
  # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
  
  idx <- which(grepl(paste0("^", comp_grp, "$"), all_subset_clean_grouped$compound_group))
  
  if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
    all_similar_compounds_idx1 <- c(all_similar_compounds_idx1, idx)
  }
  else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
    all_unique_compounds_idx1 <- c(all_unique_compounds_idx1, idx)
  }
  else {
    all_different_compounds_idx1 <- c(all_different_compounds_idx1, idx)
  }
}})

similar_compounds <- all_subset_clean_grouped[all_similar_compounds_idx1,]
other_compounds <- all_subset_clean_grouped[all_different_compounds_idx1,]
unique_compounds <- all_subset_clean_grouped[all_unique_compounds_idx1,]

# # Faster if statement
# library(data.table)
# DT <- data.table(z)
# DT[, id:=.I]
# 
# DT[, cond1:=V2!=min(V2), by=V1]

# Approach 2: Using original compound name
# all_similar_compounds_idx2 <- c()
# all_different_compounds_idx2 <- c()
# all_unique_compounds_idx2 <- c()
# 
# system.time({for (com in unique(all_subset_clean_grouped$Compound)) {
#   # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64
#   
#   # Filter 1: for the weird compound names
#   if (class(try(which(grepl(paste0("^", com, "$"), all_subset_clean_grouped$Compound)), silent = TRUE)) %in% "try-error") {
#     idx <- which(grepl(com, all_subset_clean_grouped$Compound, fixed = TRUE))
#     if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
#       # append the remove slice_df to all_subset_clean_similar
#       all_similar_compounds_idx2 <- c(all_similar_compounds_idx2, idx)
#     }
#     else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
#       all_unique_compounds_idx2 <- c(all_unique_compounds_idx2, idx)
#     }
#     else {
#       all_different_compounds_idx2 <- c(all_different_compounds_idx2, idx)
#     }
#   }
#   else {
#     if (identical(which(grepl(paste0("^", com, "$"), all_subset_clean_grouped$Compound)), integer(0))) {
#       idx <- which(grepl(com, all_subset_clean_grouped$Compound, fixed = TRUE))
#       if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
#         # append the remove slice_df to all_subset_clean_similar
#         all_similar_compounds_idx2 <- c(all_similar_compounds_idx2, idx)
#       }
#       else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
#         all_unique_compounds_idx2 <- c(all_unique_compounds_idx2, idx)
#       }
#       else {
#         all_different_compounds_idx2 <- c(all_different_compounds_idx2, idx)
#       }
#     }
#     else {
#       # Filter 2: for regular compound names
#       idx <- which(grepl(paste0("^", com, "$"), all_subset_clean_grouped$Compound))
#       if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
#         # append the remove slice_df to all_subset_clean_similar
#         all_similar_compounds_idx2 <- c(all_similar_compounds_idx2, idx)
#       }
#       else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
#         all_unique_compounds_idx2 <- c(all_unique_compounds_idx2, idx)
#       }
#       else {
#         all_different_compounds_idx2 <- c(all_different_compounds_idx2, idx)
#       }
#     }
#   }
# }})

# When include 99% of cumulative peak height, all diesel samples share 304 compounds in common
# When include 99% of cumulative peak height, all gasoline samples share 39 compounds in common
# When include 99% of cumulative peak height, all diesel composite samples share 357 compounds in common
# When include 99% of cumulative peak height, all gasoline composite samples share 248 compounds in common
# When include 99% of cumulative peak height, all IL samples share 29(22, 13) compounds in common
# After grouping compounds based on RT1, RT2, Ion1 threshold, all 31 IL samples share 13 compound "groups" in common 


# Data Summary ----------------------------------------------------------------------------------------------------
# Distribution of %Area of each compound_group --> histogram plot
ggplot(data = similar_compounds, aes(x = Percent_Area)) +
  facet_wrap(~compound_group, scales = "free_y") +
  geom_histogram(bins = 50)

# Distribution of RT1, RT2, Ion1 --> Data.table
summarydata <- similar_compounds %>%
  group_by(fuel_type, compound_group) %>%
  summarise(
    n = n(compound_group)
            )

# Distribution of compound in different fuel types
ggplot(data = similar_compounds, aes(x = compound_group)) +
  facet_wrap(~fuel_type, scales = "fixed") +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90))

# Distribution of Compound frequency in similar_compound -> Frequency bar plot



# PCA -------------------------------------------------------------------------------------------------------------
pcasubset <- similar_compounds %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound_group = factor(compound_group, levels = c(unique(compound_group)))) %>%
  # for a sample, if there are multiple occurences of a compound, then impute with mean of %Area and %Height 
  group_by(sample_name, compound_group) %>%
  # Here we collapse the duplicates compound by calculate the mean of Percent Area,
  # assuming that duplicates of similar compounds has the normal distribution
  summarise(across(Percent_Area, median)) %>% 
  # filter(sample_name %in% gas_clus2_sample) %>%
  pivot_wider(names_from = compound_group, values_from = Percent_Area) %>% # c(Percent_Area, Percent_Height)
  column_to_rownames(., var = "sample_name") # must do before input in imputePCA()

# remove columns that has less than 5 unique values, including NA as a unique value
# Aka. we remove compounds that exist in less than 3 samples ("lower bound compound filter")
remove_col <- c()
for (col in 1:ncol(pcasubset)) {
  if (length(unique(pcasubset[,col])) < 8) {
    remove_col <- c(remove_col, col)
  }
}

pcasubset_removecol <- subset(pcasubset, select = -remove_col)

# TRy imputePCA function, PCA input without NA values
PCA_impute <- imputePCA(pcasubset_removecol,
                        scale = TRUE,
                        maxiter = 2000,
                        method = "Regularized", #iterative approach-less overfitting
                        seed = 651)

# For plotting biplot later-------------------------------
subset2 <- rownames_to_column(data.frame(PCA_impute$completeObs),
                                             "sample_name")

subset_filterquantile2 <- subset2 %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 
  # %>% filter(sample_name %in% gas_clusall)

pca_input <- data.frame(PCA_impute$completeObs)
res.pca <- PCA(pcasubset, 
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca,
         addlabels = TRUE)
# Biplot
fviz_pca_biplot(res.pca,
                select.var = list(cos2 = 5),# name, # Top x active variables with the highest cos2
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = pcasubset$sample_name,
                # addEllipses=TRUE,
                dpi = 480)

# Top variables (RT1, RT2,etc.) and compounds with highest contribution
fviz_contrib(res.pca, choice = "var",
             top = 1500,
             axes = 1:2) + # contrib of var to PC1 and 2
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))

# Extract the top 1100 compounds contribute the most to PC1:PC2
var_contrib_sorted <- data.frame(var_contrib) %>%
  rownames_to_column(., var = "Compound") %>%
  mutate_at("Dim.1", funs(sort(., decreasing = TRUE))) %>% # sort descending percent_area
  mutate_at("Dim.2",funs(sort(., decreasing = TRUE)))

var_contrib_sorted <- slice_head(df, n = 1500)

# Hierarchical Clustering on Principle Components
hcpc <- HCPC(res.pca, nb.clust = -1)


# Regression Classification PCR/PLS-DA, etc. ----------------------------------------------------------------------
library(caret)
subset_filterquantile2$lab_enc <- ifelse(str_detect(subset_filterquantile2$sample_name, "DieselCom"), 1, 
                                         ifelse(str_detect(subset_filterquantile2$sample_name, "GasComp"), 2,
                                                ifelse(str_detect(subset_filterquantile2$sample_name, "D"), 3, 4)))

subset_filterquantile2 <- subset_filterquantile2 %>%
  relocate(lab_enc, .before = 2)
  
# Partitioning data set 
classification_data <- subset(subset_filterquantile2, select = -c(sample_name))
inTraining <- createDataPartition(classification_data$lab_enc,
                                  p = .60, list = FALSE)
training <- classification_data[inTraining,]
testing  <- classification_data[-inTraining,]

ctrl <- trainControl(
  method = "cv",
  number = 10,
)

model <- train(lab_enc~.,
               data = training,
               method = "pls", # "lm", pls", "lasso", rf"
               preProcess = c("center", "scale", "pca"),
               trControl = ctrl)

plot(model)

test.features <- subset(testing, select = -c(lab_enc))
test.target <- subset(testing, select = lab_enc)[,1]

predictions <- predict(model, newdata = test.features)
cor(test.target, predictions) ^ 2


var <- get_pca_var(res.pca)
var_coord <- var$coord
var_contrib <- var$contrib
# var_cos2 <- var$cos2
# corrplot(var$cos2, is.corr = FALSE)


# t-SNE clustering ------------------------------------------------------------------------------------------------
# REFERENCES VISUALIZATION: https://plotly.com/r/t-sne-and-umap-projections/
# https://distill.pub/2016/misread-tsne/
features <- subset(pcasubset, select = -c(sample_name)) 
# subset_filterquantile_similar - produced dissimilar result to PCA on the same dataset
tsne <- tsne(features,
             initial_dims = 3, 
             k = 3, 
             perplexity = 15, # Hyperparameter: perplexity < number of data points
             max_iter = 2000
             )
             # pca = FALSE, perplexity=10, theta=0.5, dims=2,
             # check_duplicates = FALSE)

pdb <- cbind(data.frame(tsne),pcasubset$sample_name)
options(warn = -1)
tsne_plot <- plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, 
               color = ~pcasubset$sample_name) %>% 
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene =list(bgcolor = "#e5ecf6"))
tsne_plot


# UMAP Clustering -------------------------------------------------------------------------------------------------
umap <- umap(features, n_components = 3, random_state = 15)

layout <- cbind(data.frame(umap[["layout"]]), pcasubset$sample_name)
umap_plot <- plot_ly(layout, x = ~X1, y = ~X2, z = ~X3, 
                color = ~pcasubset$sample_name) %>% 
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'x-axis'), 
                                   yaxis = list(title = 'y-axis'), 
                                   zaxis = list(title = 'z-axis'))) 
umap_plot

# REDUNDANT CODEs
# Selecting representative Diesel compounds for each diesel sample ---------------------------------------------------------------------------
# Get coordinates of variables
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

# Quadrant 1  - dim1 [0 to 1], dim2 [0 to -0.25] - Influencer for sample 0220F001D
sample_0220F001D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.7 & Dim.2 > -0.4 & Dim.2 < -0.15)

# Quadrant 2 - dim1 [0 to 1], dim2 [0 to 0.15] - Influencer for sample 0220F009-2D 
sample_0220F009_2D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.5 & Dim.2 > 0.05 & Dim.2 < 0.1)

# Quadrant 3 - dim1 [0 to 1], dim2 [0.15 to 0.25] - Influencer for sample 0220F009D
sample_0220F009D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0.5 & Dim.1 > 0.25 & Dim.2 > 0.3 & Dim.2 < 0.4)

# Quadrant 4  - dim1 [0 to -1], dim2 [0.6 to 0.9]- Influencer for sample 0220F005D
sample_0220F005D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.3 & Dim.1 > -0.5 & Dim.2 > 0.6 & Dim.2 < 0.9)

# Quadrant 5  - dim1 [0 to -1], dim2 [-0.5 to -0.75]- Influencer for sample 0220FDieselComp1,2
sample_0220FDieselComp1_2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.7 & Dim.2 < -0.5 & Dim.2 > -0.75)

# Quadrant 6  - dim1 [0 to -1], dim2 [-0.75 to -1]- Influencer for sample 0220FDieselComp3
sample_0220FDieselComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.15 & Dim.2 < -0.7)

# Biplot
name <- list(name = sample_0220FDieselComp3_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <- fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220FDieselComp3_influencer.png"),
              biplot_pca,
              height = 8,
              width = 15)

# Selecting representative Gasoline Composite compounds for each Gasoline sample ---------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord
sample_0220GasComp1_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.025 & Dim.1 < 0.06 & Dim.2 > 0)

sample_0220GasComp2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.85 & Dim.2 < -0.425 & Dim.2 > -0.47)

sample_0220GasComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.7 & Dim.2 < -0.5 & Dim.2 > -0.55)

# Biplot
name <- list(name = sample_0220GasComp2_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220GasComp2_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)


# Selecting representative DieselComp compounds for each dieselcomp sample ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord
# Quadrant 1  - dim1 [0 to 1], dim2 [0 to -0.25] - Influencer for sample 0220F001D
sample_0220FDieselComp1_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.75 & Dim.2 > 0.55 & Dim.2 < 0.6)

# Quadrant 2 - dim1 [0 to 1], dim2 [0 to 0.15] - Influencer for sample 0220F009-2D 
sample_0220FDieselComp2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.15 & Dim.1 > -0.19 & Dim.2 < -0.75)

# Quadrant 3 - dim1 [0 to 1], dim2 [0.15 to 0.25] - Influencer for sample 0220F009D
sample_0220FDieselComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.75 & Dim.2 > 0.35 & Dim.2 < 0.4)

# Biplot
name <- list(name = sample_0220FDieselComp3_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220FDieselComp3_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# Selecting representative Gasoline station 1, 3, 8 by gas 91 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00191_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 < -0.4 & Dim.2 > -0.55)

sample_0220F00391_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.2 < -0.4 & Dim.2 > -0.6)

sample_0220F00891_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.1 > -0.1 & Dim.2 > 0.75)

# Biplot
name <- list(name = sample_0220F00891_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00891_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)

# Selecting representative Gasoline station 1, 3, 8 by gas 87 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00187_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 > 0.15 & Dim.2 < 0.25)

sample_0220F00387_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.5 & Dim.2 <  -0.6)

sample_0220F00887_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.5 & Dim.2 > 0.66)

# Biplot
name <- list(name = sample_0220F00887_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00887_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# Selecting representative Gasoline station 1, 3, 8 by gas 89 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00189_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.15 & Dim.1 < 0.35 & Dim.2 < -0.8)

sample_0220F00389_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.2 > 0.2 & Dim.2 < 0.4)

sample_0220F00889_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 > 0.56 & Dim.2 < 0.66)

# Biplot
name <- list(name = sample_0220F00889_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00889_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# PCA --------------------------------------------------------------------
# Using Percent_Area and Percent_height value on compounds
# Diesel and DieselComp cluster
# diesel_sample <- c("0220F001D.xlsx", "0220F005D.xlsx", "0220F009-2D.xlsx", "0220F009D.xlsx",
#                    "0220FDieselComp1.xlsx", "0220FDieselComp2.xlsx", "0220FDieselComp3.xlsx")
# dieselcomp_sample <- c("0220FDieselComp1.xlsx", "0220FDieselComp2.xlsx", "0220FDieselComp3.xlsx")
# 
# # All gasoline
# gas_clusall <- c(
#   "0220F00187-2.xlsx","0220F00187-3.xlsx","0220F00187.xlsx","0220F00387.xlsx", "0220F00887.xlsx"
#   ,"0220F00189.xlsx","0220F00389.xlsx", "0220F00889.xlsx"
#   ,"0220F00191.xlsx", "0220F00391.xlsx","0220F00891.xlsx"
#   ,"0220F00894.xlsx", "0220F00587.xlsx","0220F00589.xlsx","0220F00591.xlsx",
#   "0220F00787.xlsx","0220F00789.xlsx","0220F00791.xlsx",
#   "0220F00987.xlsx", "0220F00989.xlsx", "0220F00991.xlsx"
# )
# 
# # Gasoline station 1, 3, 8 cluster - cluster 1
# gas_clus1_sample <- c(
# "0220F00187-2.xlsx","0220F00187-3.xlsx","0220F00187.xlsx","0220F00387.xlsx", "0220F00887.xlsx"
# ,"0220F00189.xlsx","0220F00389.xlsx", "0220F00889.xlsx"
# ,"0220F00191.xlsx", "0220F00391.xlsx","0220F00891.xlsx"
# ,"0220F00894.xlsx", 
# )
# 
# # Gasoline station 5, 7, 9 cluster - cluster 2
# gas_clus2_sample <- c("0220F00587.xlsx","0220F00589.xlsx","0220F00591.xlsx",
#                        "0220F00787.xlsx","0220F00789.xlsx","0220F00791.xlsx",
#                        "0220F00987.xlsx", "0220F00989.xlsx", "0220F00991.xlsx")
# 
# # Checkpoint for %Area and %Height distribution of unique compounds in all_similar_compounds
# ggplot(all_similar_compounds, aes(x = Percent_Height)) + 
#   facet_wrap(~Compound, scales = "free_y") + 
#   geom_histogram(binwidth = 0.0005) 
# Checkpoint Result: OK, it seems like both %Area and %Height have centralized distribution->it"s safe to calculate mean
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

# PCA on indi_IL_type ---------------------------------------------------------------------------------------------
# Using only Percent_height value on compounds that exist in only one sample
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

hist(all_similar_compounds$Percent_Height,
     # xlim = c(0, 0.002), 
     breaks = 3000)
quantile(all_similar_compounds$Percent_Height)
filter_quantile <- subset(all_similar_compounds, 
                          Percent_Height >  0.0005 &
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


# PCA with all_similar_compounds - Worse PC percentage explained variability than imputePCA -----------------------
subset_filterquantile_similar <- all_similar_compounds %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  # for a sample, if there are multiple occurences of a compound, then impute with mean of %Area and %Height
  group_by(sample_name, Compound) %>%
  summarise(across(Percent_Area, mean)) %>% # c(Percent_Area, Percent_Height) assuming that duplicates of similar compounds has the normal distribution
  # filter(sample_name %in% gas_clusall) %>%
  pivot_wider(names_from = Compound, values_from = Percent_Area)# c(Percent_Area, Percent_Height)

all_similar_compounds_PCA <- PCA(subset_filterquantile_similar[c(2:dim(subset_filterquantile_similar)[2])],
                                 scale.unit = TRUE,
                                 graph = FALSE)

# Scree plot
fviz_eig(subset_filterquantile_similar,
         addlabels = TRUE)

fviz_pca_var(subset_filterquantile_similar,
             col.var = "black",
             select.var = list(cos2 = 20))

fviz_pca_biplot(all_similar_compounds_PCA,
                select.var = list(cos2 = 50),# name, # list(cos2 = 140), # Top x active variables with the highest cos2
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset_filterquantile_similar$sample_name,
                # addEllipses=TRUE,
                dpi = 480)

# ggsave(paste0(getwd(), "/PCA graphs/Gasoline_station5_7_9.png"),
#        Gasoline_station5_7_9,
#        height = 8,
#        width = 15)



# Without imputePCA function, PCA input full of NA values - REALLY BAD CLUSTERING than imputePCA ------------------
# subset_filterquantile1 <- rownames_to_column(subset_filterquantile_removecol, "sample_name") 
# subset_filterquantile1$sample_name <- factor(subset_filterquantile1$sample_name, 
#                                              levels = c(unique(subset_filterquantile1$sample_name)))
# subset_filterquantilePCA <- PCA(subset_filterquantile1[c(2:dim(subset_filterquantile1)[2])], 
#                                 scale.unit = TRUE, 
#                                 graph = FALSE)
# 
# fviz_pca_biplot(subset_filterquantilePCA,
#                 select.var = list(cos2 = 50),# name, # list(cos2 = 140), # Top x active variables with the highest cos2
#                 repel = TRUE,
#                 axes = c(1,2),
#                 label = "ind",
#                 habillage = subset_filterquantile1$sample_name,
#                 # addEllipses=TRUE,
#                 dpi = 480)



# Reserve code ----------------------------------------------------------------------------------------------------
all_similar_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
colnames(all_similar_compounds) <- colnames(all_subset_clean)
all_similar_compounds$Compound <- as.character(all_similar_compounds$Compound)
all_similar_compounds$sample_name <- as.character(all_similar_compounds$sample_name)

all_different_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
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

length(unique(all_similar_compounds$Compound))
length(unique(all_unique_compounds$Compound))
# Label compound type based on chemical structure - MAYBE IRRELEVANT -----------------------------------------------------------------
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


#--------------------------------------------------------------------------------------
# 28 out of 28 gasoline samples share 25 common compounds
# more than 15 out of 28 gasoline samples share 96 common compounds
# 5 out of 5 diesel samples share 50 common compounds
# more than 15 out of 39 IL samples share 150 common compounds
# 963 compounds unique for 5 diesel compounds

# examine the cumulative peak height and peak area per sample of compounds found across 39 samples