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
# library(foreach)
# library(doSNOW)
library(writexl)
# library(rapportools)
# library(randomForest)
# library(e1071)
# library(gbm)
# library(fitdistrplus)
# library(RcppML)
# library(CAMERA)
# library(NMF)
# library(Boruta)
# library(mixOmics)
# library(VIM)
# source("RUVRand.R")

# vignette("parallel")
# options(ggrepel.max.overlaps = 300)
# set.seed(12345)
# cl <- makeCluster(8, type = "SOCK")
# registerDoSNOW(cl)

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
                 data$RT2 < (rt2 + 0.125) & data$RT2 > (rt2 - 0.125) & 
                 data$`Ion 1` < (ion1 + 0.1) & data$`Ion 1` > (ion1 - 0.1) & 
                 is.na(data$compound_group))
    if (identical(idx, integer(0))) {
      next
    }
    else {
      data[idx, "compound_group"] <- paste0("Compound_", i, ".")
      i <- i + 1
    }
    rm(rt1)
    rm(rt2)
    rm(ion1)
  }
  return(data)
}

# Filtering similar and unique compound
comp_filter <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$compound_group)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$compound_group, fixed = TRUE))
    
    if (length(unique(data[idx,]$sample_name)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$sample_name)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# Probabilistic Quotient Normalization
pqn <- function(X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  
  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }
  
  return(X.norm)
}

# RLA plots
RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
                     cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...) {
  type <- match.arg(type)
  groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
  unique.groups <- levels(groups)
  if (is.null(cols)) 
    cols <- ColList(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(inputdata))))
  for (ii in 1:length(inputdata[, 1])) 
    box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
  
  # Within groups
  if(type == "wg") {
    out_data<-data.frame()
    for (grp in unique.groups) {
      submat <- inputdata[which(inputdata[, 1] == grp), -1]
      med_vals <- apply(submat, 2, median)
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_data <- rbind(out_data, swept_mat)
    }
    # Across groups (i.e. type == "ag")
  } else  {
    med_vals <- apply(inputdata[, -1], 2, median)
    out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
  }
  
  boxplot(t(out_data),
          cex.axis=cex.axis,                 # font size
          las=las,                           # label orientation
          col=box_cols,                      # colours
          ylim=ylim,                         # y-axis range
          oma=oma,                           # outer margin size
          ...
  )
  
  abline(h=0)
}


# Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')

# Pipe operator for isolating IL types
gas_gascomp <- file_list %>% # indi_IL_file_list
  .[!str_detect(., "D")]  %>%
  .[!str_detect(., "DieselComp")] %>%
  # .[str_detect(., "GasComp")] %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")]

# Import IL samples to list
df_list <- purrr::map(gas_gascomp, read_xlsx, sheet = "Results") # indi_IL_file_list

# Filtering out column bleed and solvent --------------------------------------
df_list_clean <- purrr::map(df_list, filtering, filter_list = c("^Carbon disulfide$", 
                                                         "^Benzene$", 
                                                         "Cyclotrisiloxane..hexamethyl",
                                                         "Cyclotetrasiloxane..octamethyl",
                                                         "^Toluene$",
                                                         "^Ethylbenzene$",
                                                         "Xylene")) 



# 1st layer Normalization with TSN and Data distribution Pre-normalization ---------------------------------------------------------------------------
# data_plot_pre_norm <- list()
# for (i in 1:length(df_list_clean)) {
#   data_plot_pre_norm[[i]] <- ggplot(data = df_list_clean[[i]],
#                          aes(x = Area)) +
#     geom_histogram(bins = 100) +
#     ggtitle(indi_IL_file_list[[i]]) + 
#     scale_x_continuous(breaks=seq(0, 5000000, 1000000), limits = c(0, 5000000))
# }
# grid.arrange(grobs = data_plot_pre_norm, ncol = 5) # Data sets are all heavy left-skewed

# 1st layer Normalization with TSN
slice_df_list <- list() 

system.time({for (i in 1:length(df_list_clean)) { 
  df <- df_list_clean[[i]] %>%
    # Percentage-based normalization aka. Total Sum normalization
    mutate(Percent_Area = Area/sum(Area)) %>%
    mutate(Percent_Height = Height/sum(Height)) %>%
    arrange(desc(Percent_Height)) # Percent_Height / Percent_Area

  # # subset data based on the largest number of iteration
  # for (row_num in 1:nrow(df)) {
  #   # slice data based on condition of cumulative sum of percent_height
  #   if (sum(df[1:row_num,]$Percent_Height) > 0.99) { # Percent_Height / Percent_Area
  #     new_df <- slice_head(df, n = row_num)
  #     break
  #   }
  # }
 
  # Add sample_name column 
  slice_df_list[[i]] <- df %>% 
    mutate(sample_name = gas_gascomp[[i]]) %>% # indi_IL_file_list
    mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                              ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                     ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
}})

# Data distribution post-TSN - Data sets are all heavy left-skewed
# data_plot_TSN <- list()
# for (i in 1:length(slice_df_list)) {
#   data_plot_TSN[[i]] <- ggplot(data = slice_df_list[[i]],
#                            aes(x = Percent_Area)) +
#     geom_histogram(bins = 100) +
#     ggtitle(indi_IL_file_list[[i]])
# }
# grid.arrange(grobs = data_plot_TSN, ncol = 5)

# Grouping compounds based on RT1, RT2, Ion1 -----------------------------------------------------------------------
# Combine all subset_df together
all_data_pre_norm <- bind_rows(slice_df_list)

all_data_pre_norm_grouped <- grouping_comp(all_data_pre_norm)

# Export df for later use
# write_xlsx(all_data_pre_norm_grouped, paste0(getwd(), "/grouping_compounds_byPercent_height-240622.xlsx"))
# testing_import <- read_excel(paste0(getwd(), "/grouping_compounds.xlsx"))


# Pivot wider for other type of Data Normalization ---------------------------------------------------------------------------------
all_data_pre_norm_grouped_wider <- all_data_pre_norm_grouped %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound_group = factor(compound_group, levels = c(unique(compound_group)))) %>%
  # for a sample, if there are multiple occurences of a compound, then impute with mean of %Area and %Height 
  group_by(sample_name, compound_group) %>% # fuel_type 
  # Here we collapse the duplicates compound by calculate the mean of Percent Area,
  # assuming that duplicates of similar compounds has the normal distribution
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound_group, values_from = Percent_Area)

# View(all_data_pre_norm_grouped_wider)

# Data Normalization and data reduction [based on cumulative sum(Percent_Height)]  --------
# Median Normalization
MN <- function(data){
  #Create an empty dataframe to fill with the transformed data
  MEDIANdf <- data[,c(1:2)] 
  
  #Assign the reference profile (here arbitrarily chosen as the first profile in the dataset)
  ref <- data[1, 3:length(data)]
  append_list <- list()
  append_list[[1]] <- data[1, 3:length(data)]
  #Loop across each profile in the subsampled dataset
  for (row in 2:nrow(data)) {
    #Calculate the median ratio between the profile being normalised at row, and the reference
    ref_med <- median(na.omit(unlist(data[row, 3:length(data)]/ref)))
    #Then divide all peaks in the profile by the median ratio
    append_list[[row]] <- data[row,c(3:length(data))]/ref_med
    }
  new_mediandf <- bind_cols(MEDIANdf, bind_rows(append_list))
  #Return the normalised data
  return(new_mediandf)
}

PQN <- function(data) {
  #Create an empty data frame that will be used as the median reference
  ref <- c()
  newdata <- data.table(data[,c(2:length(data))])
  PQNdf <- all_data_pre_norm_grouped_wider[,1]
  #Then calculate the median of each peak across all sample profiles
  for (j in c(2:length(newdata))) {
    ref <- c(ref, median(newdata[[j]]))
  }
  quotients_list <- list()
  #Loop across each profile in a dataset
  for (j in c(2:length(data))) {
    quotients_list[j] <- data[,j]/ref[j]
  }
  
  #Calculate the median quotient
  median_all <- median(na.omit(unlist(quotients_list)))
  quotients_list2 <- list()
  for (i in 1:nrow(data)) {
    quotients_list2[[i]] <- data[i,]/median_all
  }
  new_PQNdf <- bind_cols(PQNdf, bind_rows(quotients_list2))
  return(new_PQNdf)
}

MN_wider_data <- MN(all_data_pre_norm_grouped_wider)

data_plot_MN <- MN_wider_data %>%
  pivot_longer(cols = c(3:length(.)),
               names_to = "compound_group",
               values_to = "Percent_Area",
               values_drop_na = TRUE
               )

ggplot(data = data_plot_MN,
       aes(x = Percent_Area)) +
  facet_wrap(~ sample_name, scales = "free_y") +
  geom_histogram(bins = 100)


coda <- function(data){
  #Create an empty dataframe to fill with the transformed data
  coda.sims <- data.frame(matrix(ncol = length(data), nrow = nrow(data)))
  #Loop across each chromatogram in the subsampled dataset
  for(i in 1:nrow(data)){
    #Build a function for calculating the geometric mean
    gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
      if(any(x < 0, na.rm = TRUE)){
        return(NaN)
      }
      if(zero.propagate){
        if(any(x == 0, na.rm = TRUE)){
          return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
      } else {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }
    }
    #calculate the sum of all peaks in a chromatogram
    g.mean <- gm_mean(data[i,c(4:length(data))])
    #Then loop across all peaks in a chromatogram
    for(j in c(4:length(data))){
      #and divide each peak by the geometric mean for that individual profile
      coda.sims[i,j] <- log(data[i,j]/g.mean)
      coda.sims[i,j][is.infinite(coda.sims[i,j])] <- 0#Then add the original individual identifiers to the normalised data
      coda.sims[,c(1:3)] <- data[,c(1:3)]
      #Assign the original column names
      colnames(coda.sims) <- names(data)
    }}
  #Finally return the normalised data
  return(coda.sims)
}
# Similar Compounds found across samples ------------------------------------------------
# Approach 1: Using compound "groups" by RT1, RT2, Ion1 Threshold
system.time({idx_list <- comp_filter(all_data_pre_norm_grouped, length(indi_IL_file_list))})


similar_compounds <- all_data_pre_norm_grouped[idx_list[[1]],][, -c(2,3,6:8)]
other_compounds <- all_data_pre_norm_grouped[idx_list[[2]],][, -c(2,3,6:8)]
unique_compounds <- all_data_pre_norm_grouped[idx_list[[3]],][, -c(2,3,6:8)]


pairwise_test <- function(df, p_val_threshold, test_choice = "t.test") {
  pairwise_test <- list()
  i <- 1
  for (com_grp in unique(df$compound_group)) {
    templist <- list()
    # iterates through every fuel_type
    for (sample1 in unique(df$sample_name)) { #  sample1 in unique(df$sample_name) / fuel1 in unique(df$fuel_type) 
      idx1 <- which(df$sample_name == sample1 & df$compound_group == com_grp) #   / df$fuel_type == fuel1
      if (length(idx1) < 2) {
        next
      }
      else {
        for (sample2 in unique(df$sample_name)) { #  sample2 in unique(df$sample_name) / fuel2 in unique(df$fuel_type)
          if (sample1 == sample2) { # sample1 == sample2 / fuel1 == fuel2
            next
          }
          else {
            idx2 <- which(df$sample_name == sample2 & df$compound_group == com_grp) #  df$sample_name == sample2 / df$fuel_type == fuel2
            if (length(idx2) < 2) {
              next
            }
            else {
                if (test_choice == "ks") {
                  templist[paste0(sample1, "-", sample2)] <- ks.test(x = df[idx1,]$Percent_Area,
                                                                 y = df[idx2,]$Percent_Area,
                                                                 alternative = "two.sided")$p.value
                }
                else if (test_choice == "mn") {
                  templist[paste0(sample1, "-", sample2)] <- wilcox.test(x = df[idx1,]$Percent_Area,
                                                                     y = df[idx2,]$Percent_Area,
                                                                     alternative = "two.sided")$p.value
                }
                else {
                  templist[paste0(sample1, "-", sample2)] <- t.test(x = df[idx1,]$Percent_Area,
                                                                y = df[idx2,]$Percent_Area)$p.value
                }
            }
          }
        }
      }
    }
    
    if (any(templist > p_val_threshold)) {
      next
    }
    else {
      pairwise_test[[i]] <- templist
      names(pairwise_test)[i] <- com_grp
    }
    i <- i + 1
  }
  return(pairwise_test)
}

# Pair wise test all data
pairwise_all_data_ttest <- pairwise_test(all_data_pre_norm_grouped, 
                                         p_val_threshold = 0.15,
                                         test_choice = "t.test")

pairwise_all_data_ks <- pairwise_test(all_data_pre_norm_grouped, 
                                         p_val_threshold = 0.15,
                                         test_choice = "ks")
pairwise_all_data_mn <- pairwise_test(all_data_pre_norm_grouped, 
                                         p_val_threshold = 0.15,
                                         test_choice = "mn")

pairwise_all_data_clean_mn <- c()
for (len in 1:length(pairwise_all_data_mn)) {
  if (length(pairwise_all_data_mn[[len]]) > 0) { # > 5
    pairwise_all_data_clean_mn <- c(pairwise_all_data_clean_mn, pairwise_all_data_mn[len])
  }
}

pairwise_all_data <- unique(names(c(pairwise_all_data_clean_ttest, 
                                    pairwise_all_data_clean_ks,
                                    pairwise_all_data_clean_mn)))

# Pair wise test similar compounds =====================================================
system.time({pairwise_similar_compounds_ttest <- pairwise_test(similar_compounds, 
                                                               p_val_threshold = 0.01,
                                                               test_choice = "t.test")})

# Combine all resulting compound groups from KS, MN and t-test
all_similar_compounds_pairwise <- unique(names(c(pairwise_similar_compounds_ks, 
                                                 pairwise_similar_compounds_mn,
                                                 pairwise_similar_compounds_ttest)))

# Filter with resulting compound groups from pair-wise comparisons and pivot wider for PCA
similar_compounds_wider <- similar_compounds %>%
  filter(., compound_group %in% all_similar_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound_group = factor(compound_group, levels = c(unique(compound_group)))) %>%
  group_by(sample_name, compound_group) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound_group, values_from = Percent_Area)

# Pair wise test Other compounds-========================================================-
# Filter only compounds that found in all 4 fuel_type
other_compounds_4_fueltype <- other_compounds %>%
  group_by(compound_group) %>%
  filter(length(unique(fuel_type)) > 1)

# BEfore imputing missing values--> pick out compound groups that most likely to has some differences between 4 fuel types
system.time({pairwise_other_compounds <- pairwise_test(other_compounds_4_fueltype,  #  
                                                       p_val_threshold = 0.01,
                                                       test_choice = "ks")})

# Remove compound that have less than 6 cross comparison
pairwise_other_compounds_clean_ks <- c()
for (len in 1:length(pairwise_other_compounds)) {
  if (length(pairwise_other_compounds[[len]]) > 5) { # > 5
    pairwise_other_compounds_clean_ks <- c(pairwise_other_compounds_clean_ks, pairwise_other_compounds[len])
  }
}

all_other_compounds_pairwise <- unique(names(c(pairwise_other_compounds_clean_ks, 
                                               pairwise_other_compounds_clean_mn,
                                               pairwise_other_compounds_clean_ttest)))

# OPTION 1: replace the missing values with LOD range ----------------------
# Create wide df contain only compounds that has some differences between 4 fuel types, but still have missing values
other_compounds_wider1 <- other_compounds_4_fueltype %>%
  filter(., compound_group %in% all_other_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound_group = factor(compound_group, levels = c(unique(compound_group)))) %>%
  group_by(sample_name, compound_group) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound_group, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

# Impute missing value with LOD
for (col in 1:ncol(other_compounds_wider1)) {
  other_compounds_wider1[which(is.na(other_compounds_wider1[,col])), col] <- runif(length(which(is.na(other_compounds_wider1[,col]))),
                                                                                   min = 0,
                                                                                   max = min(all_data_pre_norm_grouped$Percent_Area))
}

# Initiate pca input by combining imputed LOD of other compounds with filtered similar compounds
pcainput1 <- full_join(other_compounds_wider1 %>%
                         rownames_to_column(., var = "sample_name"), similar_compounds_wider) %>%
  column_to_rownames(., var = "sample_name")

# TRy PCA with these selected & imputed compound groups
res.pca1 <- PCA(pcainput1, 
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca1,
         addlabels = TRUE)
# Biplot
subset1 <- pcainput1 %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

fviz_pca_biplot(res.pca1,
                select.var = list(cos2 = 5),
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset1$sample_name,
                dpi = 900,
                title = "PCA_Biplot_LOD_Compound groups filtering via Pair-wise test at p-value 0.01")

# Hierarchical Clustering on Principle Components
hcpc <- HCPC(res.pca1, nb.clust = -1,
             metric = "euclidean",
             method = "complete",
             graph = TRUE)

# OPTION 2: REplace missing value with imputePCA ---------------------------
# Create wide df contain only compounds that has some differences between 4 fuel types, but still have missing values
other_compounds_wider2 <- other_compounds_4_fueltype %>%
  filter(., compound_group %in% all_other_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound_group = factor(compound_group, levels = c(unique(compound_group)))) %>%
  group_by(sample_name, compound_group) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound_group, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

# Impute missing values with imputePCA
imputed_other_compounds <- imputePCA(other_compounds_wider2,
                                     scale = TRUE,
                                     maxiter = 2000,
                                     method = "Regularized",
                                     seed = 123)

# COmbine compound groups selected after pairwise test of similar_compounds and other_compounds
pcainput2 <- full_join(data.frame(imputed_other_compounds$completeObs) %>%
                        rownames_to_column(., var = "sample_name"), similar_compounds_wider) %>%
  column_to_rownames(., var = "sample_name")

# TRy PCA with these selected & imputed compound groups
res.pca2 <- PCA(pcainput2, 
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca2,
         addlabels = TRUE)
# Biplot
subset2 <- pcainput2 %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

fviz_pca_biplot(res.pca2,
                select.var = list(cos2 = 5),
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset2$sample_name,
                dpi = 900,
                title = "PCA_Biplot_imputePCA_Compound groups filtering via Pair-wise test at p-value 0.01")

# Hierarchical Clustering on Principle Components
hcpc <- HCPC(res.pca2, nb.clust = -1,
             metric = "euclidean",
             method = "complete",
             graph = TRUE)

# Boxplot of data distribution of compound groups
ggplot(data = imputed_other_compounds_long 
         %>% filter(., compound_group %in% all_other_compounds_pairwise)
       , aes(fuel_type, Percent_Area, fill = fuel_type)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.25, color = "red", alpha = 0.05) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~compound_group, scales = "free_y") +
  theme()

# Data Summary after Compound Grouping ----------------------------------------------------------------------------------------------------
# Distribution of %Area of each compound_group --> histogram plot
ggplot(data = similar_compounds, aes(x = Percent_Area)) +
  facet_wrap(~compound_group, scales = "free_y") +
  geom_histogram(bins = 50)
# Distribution of RT1, RT2, Ion1
summarydata1 <- similar_compounds %>%
  group_by(fuel_type, compound_group, sample_name) %>%
  summarise(
    RT1 = as.double(unlist(across(RT1, mean))),
    RT2 = as.double(unlist(across(RT2, mean))),
    `Ion 1` = as.double(unlist(across(`Ion 1`, mean))),
    Percent_Area = as.double(unlist(across(Percent_Area, mean))) #,
    # n = n(compound_group)
            )

# Distribution of compound in different compound_group
summarydata2 <- similar_compounds %>%
  group_by(compound_group, Compound, fuel_type) %>%
  summarise(count = n(Compound)) 



# t-SNE clustering ------------------------------------------------------------------------------------------------
# REFERENCES VISUALIZATION: https://plotly.com/r/t-sne-and-umap-projections/
# https://distill.pub/2016/misread-tsne/
features <- subset(subset2, select = -c(sample_name)) # pcasubset
# subset_filterquantile_similar - produced dissimilar result to PCA on the same dataset
tsne <- tsne(features,
             initial_dims = 3, 
             k = 3, 
             perplexity = 15, # Hyperparameter: perplexity < number of data points
             max_iter = 2000
             )
             # pca = FALSE, perplexity=10, theta=0.5, dims=2,
             # check_duplicates = FALSE)

pdb <- cbind(data.frame(tsne),subset2$sample_name)
options(warn = -1)
tsne_plot <- plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, 
               color = ~subset2$sample_name) %>% 
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


# UMAP clustering -------------------------------------------------------------------------------------------------
umap <- umap(features, n_components = 3, random_state = 15)

layout <- cbind(data.frame(umap[["layout"]]), subset2$sample_name)
umap_plot <- plot_ly(layout, x = ~X1, y = ~X2, z = ~X3, 
                color = ~subset2$sample_name) %>% 
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'x-axis'), 
                                   yaxis = list(title = 'y-axis'), 
                                   zaxis = list(title = 'z-axis'))) 
umap_plot

# REDUNDANT CODEs
# Regression Classification PCR/PLS-DA, etc. ----------------------------------------------------------------------
# NEED MORE DATA TO TRAIN CLASSIFICATION
library(caret)
subset2$lab_enc <- ifelse(str_detect(subset2$sample_name, "DieselCom"), 1, 
                          ifelse(str_detect(subset2$sample_name, "GasComp"), 2,
                                 ifelse(str_detect(subset2$sample_name, "D"), 3, 4)))
# move label to the front of df
subset2 <- subset2 %>%
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




# PCA reserve code--------------------------------------------------------------------
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
# all_subset_clean / all_subset_clean_height
all_similar_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
colnames(all_similar_compounds) <- colnames(all_subset_clean)
all_similar_compounds$Compound <- as.character(all_similar_compounds$Compound)
all_similar_compounds$sample_name <- as.character(all_similar_compounds$sample_name)
all_similar_compounds$fuel_type <- as.character(all_similar_compounds$fuel_type)

all_different_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
colnames(all_different_compounds) <- colnames(all_subset_clean)
all_different_compounds$Compound <- as.character(all_different_compounds$Compound)
all_different_compounds$sample_name <- as.character(all_different_compounds$sample_name)
all_different_compounds$fuel_type <- as.character(all_different_compounds$fuel_type)

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
all_unique_compounds$fuel_type <- as.character(all_unique_compounds$fuel_type)

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


# system.time({for (comp_grp in unique(all_subset_clean_grouped$compound_group)) {
#   # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
#   
#   idx <- which(grepl(paste0("^", comp_grp, "$"), all_subset_clean_grouped$compound_group))
#   
#   if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
#     all_similar_compounds_idx1 <- c(all_similar_compounds_idx1, idx)
#   }
#   else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
#     all_unique_compounds_idx1 <- c(all_unique_compounds_idx1, idx)
#   }
#   else {
#     all_different_compounds_idx1 <- c(all_different_compounds_idx1, idx)
#   }
# }})
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
# imputed_other_compounds_long <- data.frame(imputed_other_compounds$completeObs) %>%
#   rownames_to_column(., "sample_name") %>%
#   pivot_longer(cols = c(2:length(.)),
#                names_to = "compound_group",
#                values_to = "Percent_Area") %>%
#   mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp",
#                             ifelse(str_detect(sample_name, "GasComp"), "GasComp",
#                                    ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
# 
# system.time({pairwise_imputed_other_compounds <- pairwise_test(imputed_other_compounds_long, 
#                                                                p_val_threshold = 0.1,
#                                                                test_choice = "t.test")})

# REserve code of PCA -------------------------------------------------------------------------------------------------------------
# PCA with data normalized by TSN (Percent_Area)
pcaTSN <- all_data_pre_norm_grouped_wider %>%
  column_to_rownames(., var = "sample_name") # must do before select_col and imputePCA()

# PCA with data normalized by Median Normalization
# pcaMN <- MN_wider_data %>%
#   column_to_rownames(., var = "sample_name")

# remove columns that has less than 5 unique values, including NA as a unique value
# Aka. we remove compounds that exist in less than x samples ("lower bound compound filter")
# Since Regularized approach of imputePCA drawn initial value from Gaussian distribution, we need at
# least 
# REF: https://marketing.astm.org/acton/attachment/9652/f-f77f2c0b-9bdd-43c4-b29e-a5dc68c3a4b1/1/-/-/-/-/ja17dp.pdf#:~:text=What%20is%20the%20minimum%20number%20of%20data%20points,common%20answer%20from%20most%20statistical%20professionals%20is%20%E2%80%9C30.%E2%80%9D

select_col <- c()
for (col in 1:ncol(pcaTSN)) {
  if (sum(!is.na(pcaTSN[,col])) > 30) { # select compounds that exist in every 31 samples
    select_col <- c(select_col, col)
  }
}

pcasubset_removecol <- subset(pcaTSN, select = select_col)
# selectcolcomp <- colnames(pcasubset_removecol)
# selectcolcomp %in% unique(similar_compounds$compound_group)

# Boxplot distribution
ggplot(data = pcasubset_removecol %>%
         rownames_to_column(., "sample_name") %>%
         pivot_longer(cols = c(2:length(.)),
                      names_to = "compound_group",
                      values_to = "Percent_Area") %>%
         mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                   ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                          ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))), 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound_group, scales = "free_y")

# Apply imputePCA function, since PCA input cannot have NA values
PCA_impute <- imputePCA(pcasubset_removecol,
                        scale = TRUE,
                        maxiter = 2000, # need to optimise for best max iteration
                        method = "Regularized", #iterative approach-less overfitting
                        seed = 123)

# For plotting biplot later
subset2 <- pcasubset_removecol %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

# PCA section
pca_input <- data.frame(PCA_impute$completeObs)

res.pca <- PCA(pcasubset_removecol,  #  pca_input / pcasubset
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
                habillage = subset2$sample_name, #  subset2 / pcasubset
                # addEllipses=TRUE,
                dpi = 900)

# Hierarchical Clustering on Principle Components
hcpc <- HCPC(res.pca, nb.clust = -1)

# Top variables (RT1, RT2,etc.) and compounds with highest contribution
fviz_contrib(res.pca, choice = "var",
             top = 1500,
             axes = 1:2) + # contrib of var to PC1 and 2
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))

##### Extract the top 1100 compounds contribute the most to PC1:PC2
var_contrib_sorted <- data.frame(var_contrib) %>%
  rownames_to_column(., var = "Compound") %>%
  mutate_at("Dim.1", funs(sort(., decreasing = TRUE))) %>% # sort descending percent_area
  mutate_at("Dim.2",funs(sort(., decreasing = TRUE)))

var_contrib_sorted <- slice_head(df, n = 1500)

# Iterative loop removing variables Method 1: Remove_col that have less than x number of values -----------------------------------------------------
summary_list <- list()
i <- 1

system.time({for (id in 1:length(indi_IL_file_list)) {
  
  templist <- list()
  templist <- append(templist, id)
  
  select_col <- c()
  for (col in 1:ncol(pcaTSN)) {
    if (sum(!is.na(pcaTSN[,col])) > id) { # the amount of non-NA values of compounds must > x 
      select_col <- c(select_col, col)
    }
  }
  
  pcasubset_removecol <- subset(pcaTSN, select = select_col)
  
  # Apply imputePCA function, since PCA input cannot have NA values
  PCA_impute <- imputePCA(pcasubset_removecol,
                          scale = TRUE,
                          maxiter = 2000,
                          method = "Regularized", # iterative approach-less overfitting
                          seed = 123)
  
  # PCA section
  pca_input <- data.frame(PCA_impute$completeObs)
  
  res.pca <- PCA(pca_input,
                 scale.unit = TRUE, 
                 graph = FALSE)
  
  templist <- append(templist, get_eigenvalue(res.pca)[2,3]) # get the sum of PC1 & PC2
  summary_list[[i]] <- templist
  i <- i + 1
}})

# summary_list_df <- bind_cols(summary_list)
# count(n < 7)

# Iterative loop removing variables Method 2: remove one column at a time ----------------------------------
summary_list2 <- c()
i <- 1

select_col <- c()
n <- c()
for (col in 1:ncol(pcaTSN)) {
  n <- c(n, sum(!is.na(pcaTSN[,col])))
  if (sum(!is.na(pcaTSN[,col])) > 1) { # the amount of non-NA values of compounds must > x 
    select_col <- c(select_col, col)
  }
}

pcasubset_removecol <- subset(pcaTSN, select = select_col)

# Apply imputePCA function, since PCA input cannot have NA values
PCA_impute <- imputePCA(pcasubset_removecol,
                        scale = TRUE,
                        maxiter = 2000, # need to optimise for best max iteration
                        method = "Regularized", #iterative approach-less overfitting
                        seed = 123)

# system.time({for (colnum in 1:length(pcasubset_removecol)) {
# Remove one column at a time

# Check in the number of observation of compounds in 31 sample
# n <- c()
# for (col in 4771:length(pcasubset_removecol)) {
#   n <- c(n, sum(!is.na(pcasubset_removecol[,col])))}
# min(n)
# max(n)

pca_input <- data.frame(PCA_impute$completeObs)[, -c(1:4770)] 


res.pca <- PCA(pca_input,
               scale.unit = TRUE, 
               graph = FALSE)

subset2 <- rownames_to_column(pca_input,
                              "sample_name")
subset2 <- subset2 %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name))))

fviz_pca_biplot(res.pca,
                select.var = list(cos2 = 5),
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset2$sample_name,
                dpi = 900)

hcpc <- HCPC(res.pca, nb.clust = -1)

# Boxplot of top 101 compounds values before imputePCA
# Extract compounds from pca_input 
top101compounds <- colnames(pca_input)
top101compoundsdf <- all_data_pre_norm_grouped[which(all_data_pre_norm_grouped$compound_group %in% top100compounds),]
ggplot(data = top101compoundsdf[, -c(2,3,6:8)], 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound_group, scales = "free_y")

# Boxplot of top 101 compounds values after imputePCA
ggplot(data = pca_input %>%
         rownames_to_column(., "sample_name") %>%
         pivot_longer(cols = c(2:length(.)),
                      names_to = "compound_group",
                      values_to = "Percent_Area") %>%
         mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                   ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                          ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))), 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound_group, scales = "free_y")

# summary_list2 <- c(summary_list2, get_eigenvalue(res.pca)[2,3])

# }})


# When include 99% of cumulative peak height, all diesel samples share 304 compounds in common
# When include 99% of cumulative peak height, all gasoline samples share 39 compounds in common
# When include 99% of cumulative peak height, all diesel composite samples share 357 compounds in common
# When include 99% of cumulative peak height, all gasoline composite samples share 248 compounds in common
# When include 99% of cumulative peak height, all IL samples share 29 compounds in common (b4 compound grouping)
# After grouping compounds based on RT1, RT2, Ion1 threshold, all 31 IL samples share 13 compound "groups" in common 

