install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(data.table)
library(clusterProfiler)
library(ggplot2)

# function of GO plot
# prepare log2FC file with the name: ko gene + '_' + background gene + '_FC.txt'
# for example, truA_WT_FC.txt

GO_plot_KO <- function(KO, bkg, go_univ) {
  FC_input<-paste0(KO, '_', bkg, '_FC.txt')
  
  data=fread(FC_input, header = T, data.table = F, stringsAsFactors = FALSE)
  
  all_geneid<-unique(as.character(data$'Gene'))
  
  # load univ genes (all genes)
  ctr_geneid<-all_geneid

  universe <- data
  #colnames(universe) <- c("Gene","log2FC", 'p_value', 'q_value')
  colnames(universe) <- c("Gene","log2FC", 'q_value')
  
  # filter sig genes.
  up_genes <- subset(universe, log2FC >= 1 & q_value <= 0.05)
  
  down_genes <- subset(universe, log2FC <= -1 & q_value <= 0.05)
  
  expall <- c(up_genes$Gene, down_genes$Gene)
  
  ## enrich.
  #p<0.05，Benjamini adj p，q< 0.2
  go_rich <- enricher(gene = expall, 
                      universe = ctr_geneid,
                      TERM2GENE = go_univ[c('GO', 'Gene')], 
                      TERM2NAME = go_univ[c('GO', 'Description')], 
                      pAdjustMethod = 'BH', 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  
  # dotplot
  dotplot_out<-paste0(KO, '_', bkg, '_GOBP_expFC1q0.05_showp0.05q0.2.pdf')
  pdf(dotplot_out,width=7,height=7)
  dotplot(go_rich) + # show top 10 go term: showCategory=10,font.size = 12
    labs(title = KO) + # labs(title = paste0(KO, '/WT'))
    theme(
      axis.text.x = element_text(size = 10),         # X-axis tick label font size
      axis.text.y = element_text(size = 14),         # Y-axis tick label font size
      plot.title = element_text(size = 14),          # Title font size
      legend.text = element_text(size = 10),         # Legend text font size
      legend.title = element_text(size = 10)         # Legend title font size
    )
}


## E.coli RNAseq
setwd( path to working directory )
# put 'GO_class.txt' file here.
# put 'Ec_GO_BP.txt' here. It contains all GO-annotated genes in E.coli (biological process only).
# put log2FC file here with the name: ko gene + '_' + background gene + '_FC.txt'
# for example, truA_WT_FC.txt

# load universe gid:GO and exp gid:GO.
df_gene2GO=fread('Ec_GO_BP.txt', header = T, data.table = F, stringsAsFactors = FALSE)
names(df_gene2GO) <- c('Gene', 'GO')

# load GO class file. get from：http://geneontology.org/
go_class <- read.delim('go_class.txt', header = FALSE, stringsAsFactors = FALSE)
names(go_class) <- c('GO', 'Description', 'Source')

# merge all go of genome to go_class file.
go_universe <- merge(df_gene2GO, go_class, by = 'GO', all.x = TRUE)

# load FC data.
ko_strains <- c('dusB', 'miaA', 'mnmA', 'mnmC', 'mnmG', 'rlmN', 'thiI', 'trmJ', 'truA', 'truD')

for (i in ko_strains) {
  GO_plot_KO(i, WT, go_universe)
}

ko_strains <- c('WT', 'miaA', 'mnmA', 'mnmC', 'mnmG', 'rlmN', 'thiI', 'trmJ', 'truA', 'truD')
bkg = 'WT'

for (i in ko_strains) {
  GO_plot_KO(i, bkg, go_universe)
}

#### END of E.coli
#### GO for Vibrio cholerae o1 biovar el tor str. n16961
# put 'GO_class.txt' file here.
# put 'Vc_GO_BP.txt' here. It contains all GO-annotated genes in vibrio cholerae o1 biovar el tor str. n16961 (biological process only).
# put log2FC file here with the name: ko gene + '_' + background gene + '_FC.txt'
# for example, truA_WT_FC.txt

# load universe gid:GO and exp gid:GO.
df_gene2GO=fread('Vc_GO_BP.txt', header = T, data.table = F, stringsAsFactors = FALSE)
names(df_gene2GO) <- c('Gene', 'GO')

# merge all go of genome to go_class file.
go_universe <- merge(df_gene2GO, go_class, by = 'GO', all.x = TRUE)

# load FC data.
ko_strains <- c('dusB', 'tgt', 'truA')

bkg = 'WT'

for (i in ko_strains) {
  GO_plot_KO(i, bkg, go_universe)
}

##### END of Vibrio cholerae

# retrieve genes by GO term.

library(dplyr)

## E.coli RNAseq
setwd("/Users/decrecy-lab/Desktop/tRNA speed JMB/Jill RNAseq")

# load universe gid:GO and exp gid:GO.
df_gene2GO=fread('Ec_GO_BP.txt', header = T, data.table = F, stringsAsFactors = FALSE)
names(df_gene2GO) <- c('Gene', 'GO')

# load GO class file. get from：http://geneontology.org/
go_class <- read.delim('go_class.txt', header = FALSE, stringsAsFactors = FALSE)
names(go_class) <- c('GO', 'Description', 'Source')

# merge all go of genome to go_class file.
go_universe <- merge(df_gene2GO, go_class, by = 'GO', all.x = TRUE)

KO = 'truA'
bkg = 'WT'

FC_input<-paste0(KO, '_', bkg, '_FC.txt')

data=fread(FC_input, header = T, data.table = F, stringsAsFactors = FALSE)

all_geneid<-unique(as.character(data$'Gene'))

# load univ genes (all genes)
ctr_geneid<-all_geneid
universe <- data
#colnames(universe) <- c("Gene","log2FC", 'p_value', 'q_value')
colnames(universe) <- c("Gene","log2FC", 'q_value')

# filter sig genes.
up_genes <- subset(universe, log2FC >= 1 & q_value <= 0.05)
down_genes <- subset(universe, log2FC <= -1 & q_value <= 0.05)

expall <- c(up_genes$Gene, down_genes$Gene)

## enrich.
#p<0.05，Benjamini adj p，q< 0.2
go_rich <- enricher(gene = expall, 
                    universe = ctr_geneid,
                    TERM2GENE = go_universe[c('GO', 'Gene')], 
                    TERM2NAME = go_universe[c('GO', 'Description')], 
                    pAdjustMethod = 'BH', 
                    pvalueCutoff = 0.05 ,
                    qvalueCutoff = 0.2)

# Convert results to a data frame
go_results <- as.data.frame(go_rich)

# Look at the columns to identify where GO terms and gene lists are stored
#head(go_results)
write.table(go_results, sep = "#")


#### retrive gene names by a specific GO term.

# by the order of IDs in the result.
go_term <- "GO:0071973"

genes_of_interest <- go_results %>%
  filter(ID == go_term) %>%
  select(geneID)

# Retrieve the gene list
gene_list <- strsplit(genes_of_interest$geneID, "/")[[1]]

# read the ID mapping file.
idmap = fread('Ec_IDmapping.txt', header = T, data.table = F, stringsAsFactors = FALSE)

# Print the gene Blatter list.
print(gene_list)

# Print corresponding gene names.
matching_indices <- match(gene_list, idmap$Blattner)

# Extract the values in column 2 corresponding to the matched rows
genenames <- idmap$Gene[matching_indices]

# Print the result, maintaining the order of keys in key_list
print(genenames)


