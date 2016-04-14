# https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression#

# 1 Mouse ID 
# 2..78 Values of expression levels of 77 proteins; the names of proteins are followed by â€œ_nâ€ indicating that they were measured in the nuclear fraction. For example: DYRK1A_n 
# 79 Genotype: control (c) or trisomy (t) 
# 80 Treatment type: memantine (m) or saline (s) 
# 81 Behavior: context-shock (CS) or shock-context (SC) 
# 82 Class: c-CS-s, c-CS-m, c-SC-s, c-SC-m, t-CS-s, t-CS-m, t-SC-s, t-SC-m 

library(dplyr);
library(gplots);
library(RColorBrewer);


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


# STEP 2: Initialize the heatmap arrangement
# ==========================================

# Define the color range for the actual expression levels
my_palette <- colorRampPalette(c("green", "yellow", "orange", "red"))(n = 399);


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should change the ranges to better fit the scaling you wish
# Currently ranges are:
#  (0.0 -  0.49)  -> green
#  (0.5 -  0.79)  -> yellow
#  (0.8 -  0.99)  -> orange
#  (1.0 - 10.49)  -> red
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define the color breaks that would better represent the different values
col_breaks = c(seq(0  ,  0.49,   length=100),    # for green
               seq(0.5,  0.79,   length=100),    # for yellow
               seq(0.8,  0.99,   length=100),    # for orange
               seq(1  , 10,      length=100));   # for red


# STEP 3: Define the clustering methods
# ==========================================

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select both the distance and the clustering method
# used for better results
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clustering or rows (i.e. tissues/samples)
# Distance choices: euclidean (default), maximum, canberra, binary, minkowski, manhattan
# Cluster method choices: complete (default), single, average, mcquitty, median, centroid, ward.D, ward.D2
row_distance = dist(as.matrix(MouseDataClean)[c(1:nrow(MouseDataClean)),], method = "euclidean");
row_cluster = hclust(row_distance, method = "complete");


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# You should select both the distance and the clustering method
# used for better results
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clustering of columns (i.e. protein expressions)
# Distance choices: euclidean (default), maximum, canberra, binary, minkowski, manhattan
# Cluster method choices: complete (default), single, average, mcquitty, median, centroid, ward.D, ward.D2
col_distance = dist(t(as.matrix(MouseDataClean)[c(1:nrow(MouseDataClean)),]), method = "euclidean");
col_cluster = hclust(col_distance, method = "complete");



# STEP 4: Create heatmap
# ==========================================

# Creates the png image file
png("MouseProteinExpressionHeatmap.png",   # create PNG for the heat map        
    width  = 8000,
    height = 6000,
    res    = 300,                               # 300 pixels per inch
    pointsize = 8);                             # smaller font size

# Define a color for each of the 8 different cases
colorList = c("gray", "blue", "lightsalmon", "orchid", "skyblue", "black", "green", "chartreuse4", "burlywood");
rowCategories <- c(rep(colorList[1], nrow(MouseData_cCSs)),   # c-CS-s
                   rep(colorList[2], nrow(MouseData_cCSm)),   # c-CS-m
                   rep(colorList[3], nrow(MouseData_cSCs)),   # c-SC-s
                   rep(colorList[4], nrow(MouseData_cSCm)),   # c-SC-m
                   rep(colorList[5], nrow(MouseData_tCSs)),   # t-CS-s
                   rep(colorList[6], nrow(MouseData_tCSm)),   # t-CS-m
                   rep(colorList[7], nrow(MouseData_tSCs)),   # t-SC-s
                   rep(colorList[8], nrow(MouseData_tSCm))    # t-SC-m
);
classNames <- c("c-CS-s", "c-CS-m", "c-SC-s", "c-SC-m", "t-CS-s", "t-CS-m", "t-SC-s", "t-SC-m");

# Constructs the heatmap with the selected parameters
heatmap.2(as.matrix(MouseDataClean)[c(1:nrow(MouseDataClean)),],                    # data for heatmap
          cellnote = as.matrix(MouseDataClean)[c(1:nrow(MouseDataClean)),],         # same data set for cell labels
          main = "Mouse Protein Expression", # heat map title
          notecol="black",                   # change font color of cell labels to black
          density.info="histogram",          # turns off density plot inside color legend ("histogram","density","none"),
          trace="none",                      # turns off trace lines inside the heat map (row, column, both, none)
          tracecol="cyan",                   # character string giving the color for "trace" line
          margins =c(50,50),                 # widens margins around plot
          col=my_palette,                    # use on color palette defined earlier 
          dendrogram="both",                 # character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms
          Rowv = as.dendrogram(row_cluster), # apply default clustering method
          Colv = as.dendrogram(col_cluster), # apply default clustering method
          RowSideColors = rowCategories,     # grouping row-variables into different categories
          breaks=col_breaks                 # enable color transition at specified limits
);

# Include legend in figure
legend("topright",                                         # location of the legend on the heatmap plot
       legend = classNames,                                # category labels
       col = colorList[1:length(unique(rowCategories))],  # color key
       lty= 1,                                             # line style
       lwd = 10                                            # line width
);

# Close the file
dev.off()
