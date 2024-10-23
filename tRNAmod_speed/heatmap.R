setwd("C:/Users/yuanyifeng/Desktop/tRNA speed JMB/RNAseq Jill")

#install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#install.packages("circlize")

library(ComplexHeatmap)
library(circlize)
library(grid)

file = "heatmap_select.txt"

data <- read.table(file, header = TRUE, sep = "\t")

# use the 1st column as row name.
rownames(data) <- data$gene

# 2. Remove the first column (since it is now used as row names).
data <- data[, -1]

# 3. Convert the remaining columns to numeric. It will lose the row names.
data_numeric <- as.data.frame(lapply(data, as.numeric))

# keep the row names.
rownames(data_numeric) <- rownames(data)

# convert table to matrix.
data_matrix <- as.matrix(data_numeric)

# 4. for row colors.
# Define different lists of row names.
list_chemotaxis <- c("aer", "cheA", "cheB", "cheR", "cheW", "cheY", "cheZ", "tap", 
                     "tar", "trg")
list_motility <- c("flgA", "flgB", "flgC", "flgD", "flgE", "flgF", "flgG", "flgH", "flgI", "flgJ", 
           "flgK", "flgL", "flgM", "flgN", "flhA", "flhB", "flhC", "flhD", "flhE", "fliA", 
           "fliD", "fliE", "fliF", "fliG", "fliH", "fliI", "fliJ", "fliK", "fliL", "fliM", 
           "fliN", "fliO", "fliP", "fliQ", "fliR", "fliS", "fliT", "motA", "motB", "pdeH", "ycgR")
list_lipid <- c("dhaK", "dhaL", "dhaM")
list_stress <- c("allR", "entS", "yghG", "yhaC", "yiaA", "yohC", "katE", "dps", "elaB", "evgS", "phoU")
list_pion <- c("pstA", "pstB", "pstC", "pstS")
list_his <- c("hisA", "hisB", "hisC", "hisD", "hisF", "hisG", "hisH")
list_iron <- c("tonB", "cirA", "exbB", "exbD", "fepA", "fepC", "fepD", "fepG", "fes", "fhuA", 
                "fhuB", "fhuC", "fhuD", "fiu")
list_aa <- c("dsdA", "alaE", "cyuP", "dsdX", "glnH", "glnP", "glnQ", "osmF", "yehY")
list_met <- c("metB", "metE", "metF", "metR")
list_carb <- c("gpmM", "malP", "malQ", "nanA", "agaI", "pflB", "uidA", "ycjM", "ygaQ", 
               "cmtA", "malX", "gatA", "gatB", "gatC", "gatZ", "galE", "mglA", "mglB", 
               "mglC", "ptsG", "glpK", "glpQ", "glpT", "glpX", "fbaB", "glk", "pfkA", 
               "manX", "manY", "manZ", "aldB", "glpD", "srlA", "srlB", "srlE", "amyA")
list_defense <- c("casA", "casB", "casC", "casD", "casE")

# Create a vector for row annotations based on the lists
# Assign background colors using hex codes based on the lists
row_colors <- ifelse(rownames(data_matrix) %in% list_chemotaxis, "#1854FC",
                ifelse(rownames(data_matrix) %in% list_lipid, "#83FDEE",
                  ifelse(rownames(data_matrix) %in% list_motility, "#47D359",
                    ifelse(rownames(data_matrix) %in% list_stress, "#B2D1F4",
                      ifelse(rownames(data_matrix) %in% list_pion, "#CC99FF",
                        ifelse(rownames(data_matrix) %in% list_his, "#CF63CA",
                          ifelse(rownames(data_matrix) %in% list_iron, "#E9A9AE", 
                            ifelse(rownames(data_matrix) %in% list_aa, "#F2ED8E", 
                              ifelse(rownames(data_matrix) %in% list_met, "#BDFD9D",
                                ifelse(rownames(data_matrix) %in% list_carb, "#F74343", 
                                  ifelse(rownames(data_matrix) %in% list_defense, "#FFC000", "white")))))))))))


# Create row annotation for the background colors (custom color function)
row_annot <- rowAnnotation(df = data.frame(Group = factor(rownames(data_matrix))),
                           col = list(Group = setNames(row_colors, rownames(data_matrix))),
                           show_legend = FALSE# Hide the annotation legend
                           )

# Define the output TIFF file with custom dimensions (e.g., 8x6 inches)
tiff("heatmap_output.tiff", width = 6, height = 8, units = "in", res = 300)

Heatmap(data_matrix, 
        name = "log2 FoldChange", 
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),  # Color scale
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        show_row_dend = FALSE,  # Hide the row dendrogram
        column_title = "Strains", 
        row_title = "Genes",
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean",
        row_names_gp = gpar(fontsize = 4, fontface = "italic"),  # Change row label font size to 10
        column_names_gp = gpar(fontface = "bold.italic"),  # Apply bold italic to column labels
        right_annotation = row_annot  # Add row annotation to the heatmap
        )

# Close the TIFF device to save the plot
dev.off()

