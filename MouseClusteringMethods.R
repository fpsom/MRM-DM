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
names(MouseDataRaw);
summary(MouseDataRaw);


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
  select(-MouseID, -Genotype, -Treatment, -Behavior);


MouseData_cCSm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-CS-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_cSCs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-SC-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_cSCm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "c-SC-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_tCSs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-CS-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_tCSm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-CS-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_tSCs <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-SC-s") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);

MouseData_tSCm <- MouseDataRaw %>%
  na.omit() %>%
  filter(class == "t-SC-m") %>%
  select(-MouseID, -Genotype, -Treatment, -Behavior);


# Finally, join all subsets in the final set. The label attribute is assigned as the
# row name, and dropped as an independent attribute
# Final construct: each row has 77 attributes (for the 77 proteins)
MouseDataClean <- bind_rows(MouseData_cCSs, MouseData_cCSm, MouseData_cSCs, MouseData_cSCm, MouseData_tCSs, MouseData_tCSm, MouseData_tSCs, MouseData_tSCm);
rownames(MouseDataClean) <- MouseDataClean$rowNamesInfo; 
MouseData <- select(MouseDataClean, -rowNamesInfo, -class);

head(MouseData);
summary(MouseData);

# STEP 2: Perform k-means clustering
# ===================================


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select:
#     a. number of centers, and
#     b. specific attributes
# that will better facilitate clustering of instances
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

kmeansK = 8;
clusters <- kmeans(MouseData, centers = kmeansK, nstart = 15, algorithm = "Hartigan-Wong");


png("KmeansMouseProteinExpression.png", # create PNG for the plot        
    width     = 8000,                   # set the width of the image in pixels
    height    = 6000,                   # set the height of the image in pixels
    res       = 300,                    # set the resolutions to 300 pixels per inch
    pointsize = 18);                     # set the size of any letters/text

# Plots the attributes DYRK1A_N and BDNF_N (columns 1 and 3 respectively) using the cluster id as different colors
plot(MouseData$DYRK1A_N, MouseData$BDNF_N, col = clusters$cluster);

# Identifies the centers for the selected columns (attributes), i.e. 1 and 3
centers <-subset(clusters$centers, select = c(1,3));

# Prints out the centroids as crosses on the figure
points(centers, col = 1:kmeansK, pch= 3, cex= 3, lwd= 3);

legend("topright",            # location of the legend on the heatmap plot
       legend = 1:kmeansK,    # category labels
       col = 1:kmeansK,       # 
       lty= 1,                # line style
       lwd = 10               # line width
);

# Closes the devices (i.e. the "print to PNG" process)
dev.off()


# Evaluate clustering  using SSE
cat("\nSEE Results:\n===\n\n")
cat(paste("  SSE Between Clusters  : ", clusters$betweenss, "\n",sep=" "));
cat(paste("SSE Total within-cluster: ", clusters$tot.withinss, "\n",sep=" "));
cat(paste("      SSE Total         : ", clusters$totss, "\n",sep=" "));

table(MouseDataClean$class, clusters$cluster)




# STEP 3: Perform hierarchical clustering
# ========================================

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select the clustering method
# that will better facilitate clustering of instances
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Construct the dissimilarity structure for the dataset
distanceMouseData <- dist(MouseData, method = "euclidean");

# perform the hierarchical clustering with "average" as the method
hc<-hclust(distanceMouseData, method="average");


library(dendextend)

png("HierarchicalClusteringMouseProteinExpression.png",     # create PNG for the plot        
    width     = 8000,                                       # set the width of the image in pixels
    height    = 6000,                                       # set the height of the image in pixels
    res       = 300,                                        # set the resolutions to 300 pixels per inch
    pointsize = 5);                                         # set the size of any letters/text

# Convert the hierarchical clustering to a dendrogram (just for visualisation)
dend <- as.dendrogram(hc)

# Create the coding for each class
groupCodes <- c(rep("cCSs", nrow(MouseData_cCSs)), rep("cCSm", nrow(MouseData_cCSm)), rep("cSCs", nrow(MouseData_cSCs)), rep("cSCm", nrow(MouseData_cSCm)),
                rep("tCSs", nrow(MouseData_tCSs)), rep("tCSm", nrow(MouseData_tCSm)), rep("tSCs", nrow(MouseData_tSCs)), rep("tSCm", nrow(MouseData_tSCm)));

# Assign different colors to different classes
colorCodes <- c(cCSs="red", cCSm="orange", cSCs="yellow", cSCm="purple",
                tCSs="blue", tCSm="darkgreen", tSCs="darkgrey", tSCm="green");

# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

# Plotting the new dendrogram
plot(dend)

legend("topright",                                                                    # location of the legend on the plot
       legend = c("cCSs", "cCSm", "cSCs", "cSCm", "tCSs", "tCSm", "tSCs", "tSCm"),    # category labels
       col = colorCodes,                                                              # color codes
       lty= 1,                                                                        # line style
       lwd = 10                                                                       # line width
);

dev.off()


