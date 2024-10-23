install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(data.table)
library(clusterProfiler)
library(ggplot2)
setwd("/Users/decrecy-lab/Desktop/tRNA speed JMB/Jill RNAseq")

GO_plot_KO <- function(KO, go_univ) {
  FC_input<-paste0(KO, '_WT_FC.txt')
  
  data=fread(FC_input, header = T, data.table = F, stringsAsFactors = FALSE)
  
  all_geneid<-unique(as.character(data$'Gene'))
  
  # load univ genes (all genes)
  ctr_geneid<-all_geneid

  universe <- data
  colnames(universe) <- c("Gene","FC", 'p_value', 'q_value')
  
  # filter sig genes.
  up_genes <- subset(universe, FC >= 1 & p_value <= 0.05)
  
  down_genes <- subset(universe, FC <= -1 & p_value <= 0.05)
  
  expall <- c(up_genes$Gene, down_genes$Gene)
  
  ## enrich.
  #p<0.05，Benjamini adj p，q< 0.2
  go_rich <- enricher(gene = expall, 
                      universe = ctr_geneid,
                      TERM2GENE = go_univ[c('GO', 'Gene')], 
                      TERM2NAME = go_univ[c('GO', 'Description')], 
                      pAdjustMethod = 'BH', 
                      pvalueCutoff = 0.05 ,
                      qvalueCutoff = 0.2)
  
  # dotplot
  dotplot_out<-paste0(KO, '_WT_GOBP_expFC1p0.05.pdf')
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

## 1. load universe gid:GO and exp gid:GO.
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
  GO_plot_KO(i, go_universe)
}

GO_plot_KO('dusB', go_universe) ; dev.off()

GO_plot_KO('miaA', go_universe) ; dev.off()

GO_plot_KO('mnmA', go_universe) ; dev.off()

GO_plot_KO('mnmC', go_universe) ; dev.off()

GO_plot_KO('mnmG', go_universe) ; dev.off()

GO_plot_KO('rlmN', go_universe) ; dev.off()

GO_plot_KO('thiI', go_universe) ; dev.off()

GO_plot_KO('trmJ', go_universe) ; dev.off()

GO_plot_KO('truA', go_universe) ; dev.off()

GO_plot_KO('truD', go_universe) ; dev.off()





###################

KO='truA'
FC_input<-paste0(KO, '_WT_FC.txt')

data=fread(FC_input, header = T, data.table = F, stringsAsFactors = FALSE)

all_geneid<-unique(as.character(data$'Gene'))
# load ctrl mod genes (all genes)
ctr_geneid<-all_geneid
length(ctr_geneid)


#go_local(mod)

universe <- data
colnames(universe) <- c("Gene","FC", 'p_value', 'q_value')

# filter sig genes
#up_genes <- subuniverse[universe$FC >= 1, ]

#down_genes <- universe[universe$FC <= -1, ]

up_genes <- subset(universe, FC >= 1 & p_value <= 0.05)

down_genes <- subset(universe, FC <= -1 & p_value <= 0.05)

exp <- c(up_genes$Gene, down_genes$Gene)


#exp_up <- up_genes$Gene
#exp_down <- down_genes$Gene

## enrich.
#p<0.05，Benjamini adj p，q< 0.2
go_rich <- enricher(gene = exp , 
                    universe = ctr_geneid,
                    TERM2GENE = go_universe[c('GO', 'Gene')], 
                    TERM2NAME = go_universe[c('GO', 'Description')], 
                    pAdjustMethod = 'BH', 
                    pvalueCutoff = 0.05 ,
                    qvalueCutoff = 0.2)

# 1. dotplot
dotplot_out<-paste0(KO, '_WT_GOBP_expFC1p0.05.pdf')
pdf(dotplot_out,width=7,height=7)
dotplot(go_rich) + # show top 10 go term: showCategory=10,font.size = 12
  labs(title = KO) +
  theme(
    axis.text.x = element_text(size = 10),         # X-axis tick label font size
    axis.text.y = element_text(size = 14),         # Y-axis tick label font size
    plot.title = element_text(size = 10),          # Title font size
    legend.text = element_text(size = 10),         # Legend text font size
    legend.title = element_text(size = 10)         # Legend title font size
  )
dev.off()

# # 2. enrichment GO plot
# ## Add similarity matrix to the termsim slot of enrichment result
# go_rich <- enrichplot::pairwise_termsim(go_rich)
# ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# clusterplot_out<-paste0(mod, '_up_GOclusterplot.pdf')
# pdf(clusterplot_out,width=7,height=5)
# emapplot(go_rich, showCategory=10)
# dev.off()
