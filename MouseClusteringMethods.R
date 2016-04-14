# https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression#

# 1 Mouse ID 
# 2..78 Values of expression levels of 77 proteins; the names of proteins are followed by â€œ_nâ€ indicating that they were measured in the nuclear fraction. For example: DYRK1A_n 
# 79 Genotype: control (c) or trisomy (t) 
# 80 Treatment type: memantine (m) or saline (s) 
# 81 Behavior: context-shock (CS) or shock-context (SC) 
# 82 Class: c-CS-s, c-CS-m, c-SC-s, c-SC-m, t-CS-s, t-CS-m, t-SC-s, t-SC-m 

library(dplyr);
library(gplots);


MouseDataRaw <- read.csv(file="data/Data_Cortex_Nuclear.csv",head=TRUE,sep=";");

summary(MouseDataRaw);
names(MouseDataRaw);
attributes(MouseDataRaw);


# STEP 1: Restructure data for the clustering process
# ====================================================
# First, include an extra column that will correspond to the label of each instance
# The label is constructed through concantenation of two existing fields:
# MouseID and class  (e.g.:   )
MouseDataRaw$rowNamesInfo <- paste(MouseDataRaw$MouseID, MouseDataRaw$class, sep="   ");

# Second, split the dataset for the different classes (i.e. 8 datasets)
# Each subset contains 78 attributes: the 77 protein expression levels, and the label
MouseData_cCSs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-CS-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);


MouseData_cCSm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-CS-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_cSCs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-SC-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_cSCm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-SC-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_tCSs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-CS-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_tCSm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-CS-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_tSCs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-SC-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);

MouseData_tSCm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-SC-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior, -class);


# Finally, join all subsets in the final set. The label attribute is assigned as the
# row name, and dropped as an independent attribute
# Final construct: each row has 77 attributes (for the 77 proteins)
MouseDataClean <- bind_rows(MouseData_cCSs, MouseData_cCSm, MouseData_cSCs, MouseData_cSCm, MouseData_tCSs, MouseData_tCSm, MouseData_tSCs, MouseData_tSCm);
rownames(MouseDataClean) <- MouseDataClean$rowNamesInfo; 
MouseDataClean <- select(MouseDataClean, -rowNamesInfo);



# STEP 2: Perform k-means clustering
# ===================================


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select:
#     a. number of centers, and
#     b. specific attributes
# that will better facilitate clustering of instances
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Creates the png image file
png("KmeansMouseProteinExpression.png",          # create PNG for the heat map        
    width     = 8000,
    height    = 6000,
    res       = 300,                             # 300 pixels per inch
    pointsize = 8)                               # smaller font size

# Performs k-mean clustering using k = 16
clusters <- kmeans(MouseDataClean, centers = 16);

# Plots the attributes DYRK1A_N and BDNF_N (columns 1 and 3 respectively)
# using the cluster id as different colors
plot(MouseDataClean$DYRK1A_N, MouseDataClean$BDNF_N, col = clusters$cluster);

# Identifies the centers for the selected columns (attributes), i.e. 1 and 3
centers <-subset(clusters$centers, select = c(1,3));
# Prints out the centroids as crosses on the figure
points(centers, col = 1:8, pch= 3, cex= 3, lwd= 3);
dev.off()


# STEP 3: Perform hierarchical clustering
# ========================================

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select the clustering method
# that will better facilitate clustering of instances
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Creates the png image file
png("HierarchicalClusteringMouseProteinExpression.png",          # create PNG for the heat map        
    width     = 8000,
    height    = 6000,
    res       = 300,                             # 300 pixels per inch
    pointsize = 8)                               # smaller font size

# perform the hierarchical clustering with "ave" (i.e. average) as the method
hc<-hclust(dist(MouseDataClean), method="ave");
plot(hc);

dev.off()


