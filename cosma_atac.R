##packages and wd
library(PISCES)
setwd('C://Users/lvlah/linux/ac_lab/cosma/')
source('cosma_functions.R')
library(pheatmap)
library(RColorBrewer)
set.seed(343)

## load gene sets
###############
## davide sets 1
klf16.promoter <- as.list(read.table('atac/BAZ2B.KLF16.promoter.genes', stringsAsFactors = FALSE))
klf5.promoter <- as.list(read.table('atac/BAZ2B.KLF5.promoter.genes', stringsAsFactors = FALSE))
sp1.promoter <- as.list(read.table('atac/BAZ2B.SP1.promoter.genes', stringsAsFactors = FALSE))
sp2.promoter <- as.list(read.table('atac/BAZ2B.SP2.promoter.genes', stringsAsFactors = FALSE))
sp3.promoter <- as.list(read.table('atac/BAZ2B.SP3.promoter.genes', stringsAsFactors = FALSE))
## davide sets 2
shared.motifs <- read.table('atac/BAZ2B_shared_vs_BAZ2B+PROG_unique_clusters_motif_names.tab', stringsAsFactors = FALSE, sep = '\t')
unique.motifs <- read.table('atac/BAZ2B_unique_vs_all_PROG_notin_LUF_clusters_motif_names.tab', stringsAsFactors = FALSE, sep = '\t')
## read atac short list
tf.mat <- read.csv('atac/atac-tf-short-list.csv', stringsAsFactors = FALSE)
tf.list <- unlist(tf.mat, use.names = FALSE)
tf.list <- tf.list[which(tf.list != '')]
tf.set <- GeneListConvert(tf.list, species = 'human', 'gn', 'ensg')
###############

## convert
## load network and matrices
k10rn.net <- readRDS('data/k10rnv_pruned.rds')
vip.mat <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
cpm.mat <- readRDS('data/exp-samps/exp-pops_cpm.rds')
## load classificaiton objects
model.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-mvc.rds')
model.pm <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-mvc.rds')
## sample filtration
cell.source <- substr(colnames(vip.mat), 1, 5)
ak11.samps <- colnames(vip.mat)[which(cell.source == 'AK011')]
ak12.samps <- colnames(vip.mat)[which(cell.source == 'AK012')]
vip.mat <- vip.mat[, c(ak11.samps, ak12.samps)]
cpm.mat <- cpm.mat[, c(ak11.samps, ak12.samps)]

## build mr sets
###############
## ak12 v ak11 mrs
class.vect <- substr(colnames(vip.mat), 1, 5); names(class.vect) <- colnames(vip.mat)
ak12v11.mrs <- BTTestMRs(vip.mat, class.vect)
###############

## check network enrichment
###############
mr.set <- names(ak12v11.mrs$AK012[1:50])
gene.set <- tf.set
file.name <- 'tf-short-list_12v11'
## modify gene set
gene.set <- unlist(gene.set, use.names = FALSE)
gene.set <- gsub("\\..*", "", gene.set)
gene.set <- intersect(gene.set, rownames(cpm.mat))
## get targets
mr.targets <- lapply(k10rn.net[mr.set], function(x) { names(x$tfmode) })
mr.targets <- unlist(mr.targets, use.names = FALSE)
mr.targets <- unique(mr.targets)
## build target vect
target.vect <- rep('non.target', nrow(cpm.mat)); names(target.vect) <- rownames(cpm.mat)
target.vect[mr.targets] <- 'target'
## build atac vect
atac.vect <- rep('non.atac', nrow(cpm.mat)); names(atac.vect) <- rownames(cpm.mat)
atac.vect[gene.set] <- 'atac'
## build table, perform fischer test
table.obj <- table(target.vect, atac.vect)
f.pVal <- fisher.test(table.obj, alternative = 'less')$p.value; print(f.pVal)
## save
save.obj <- list('table' = table.obj, 'p.val' = f.pVal)
saveRDS(save.obj, file = paste('atac/enrichment-tables/', file.name, '_target-enrichment.rds', sep = ''))
###############

## get protein activity + p-value (shared)
###############
## modify to list + convert
shared.list <- shared.motifs[,2]
shared.list <- unlist(lapply(shared.list, function(x) { unlist(strsplit(x, ',')) }))
shared.list <- unlist(lapply(shared.list, function(x) { unlist(strsplit(x, '::')) }))
shared.list <- gsub("\\(.*", "", shared.list)
shared.ensg <- GeneListConvert(shared.list, 'human', 'gn', 'ensg')
## get p-values 
p.vals <- apply(vip.mat[intersect(shared.ensg, rownames(vip.mat)),], 1, function(x) { 
  wilcox.test(x[ak12.samps], x[ak11.samps], alternative = 'greater')$p.value})
avg.vip <- rowMeans(vip.mat[intersect(shared.ensg, rownames(vip.mat)), ak12.samps])
## make table and save
shared.table <- data.frame('avg.vip' = avg.vip, 'p.val' = p.vals)
write.table(shared.table, file = 'atac/shared-motifs_vip.tsv', sep = '\t', row.names = TRUE, col.names = TRUE,
            quote = FALSE)
### HEATMAP
# cell order
annot.df <- data.frame('Class' = model.mvc[c(ak11.samps, ak12.samps)],
                       'Population' = substr(c(ak11.samps, ak12.samps), 1, 5))
annot.df <- annot.df[order(annot.df$Population, annot.df$Class),]; cell.order <- rownames(annot.df)
# colors
pact.col <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))
class.colors <- ClusterColors(4); names(class.colors) <- sort(unique(model.mvc))
pop.colors <- c('lightgrey', 'black'); names(pop.colors) <- c('AK011', 'AK012')
annot.color <- list('Class' = class.colors, 'Population' = pop.colors)
# plot.dat / breaks
plot.dat <- vip.mat[intersect(shared.ensg, rownames(vip.mat)), cell.order]
plot.dat <- GeneNameConvert(plot.dat, 'human', 'ensg', 'gn')
mat.breaks <- QuantileBreaks(plot.dat, 100)
# PLOT
dev.off()
jpeg('atac/shared-motifs_vip-heatmap.jpg', width = 1200, height = 1200)
pheatmap(plot.dat, main = 'Shared Motif Activity', fontsize = 20,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE,  show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = pact.col(length(mat.breaks) - 1))
dev.off()
###############

## get protein activity + p-value (unique)
###############
## modify to list + convert
unique.list <- unique.motifs[,2]
unique.list <- unlist(lapply(unique.list, function(x) { unlist(strsplit(x, ',')) }))
unique.list <- unlist(lapply(unique.list, function(x) { unlist(strsplit(x, '::')) }))
unique.list <- gsub("\\(.*", "", unique.list)
unique.ensg <- GeneListConvert(unique.list, 'human', 'gn', 'ensg')
## get p-values 
p.vals <- apply(vip.mat[intersect(unique.ensg, rownames(vip.mat)),], 1, function(x) { 
  wilcox.test(x[ak12.samps], x[ak11.samps], alternative = 'greater')$p.value})
avg.vip <- rowMeans(vip.mat[intersect(unique.ensg, rownames(vip.mat)), ak12.samps])
## make table and save
unique.table <- data.frame('avg.vip' = avg.vip, 'p.val' = p.vals)
write.table(unique.table, file = 'atac/unique-motifs_vip.tsv', sep = '\t', row.names = TRUE, col.names = TRUE,
            quote = FALSE)
### HEATMAP
# cell order
annot.df <- data.frame('Class' = model.mvc[c(ak11.samps, ak12.samps)],
                       'Population' = substr(c(ak11.samps, ak12.samps), 1, 5))
annot.df <- annot.df[order(annot.df$Population, annot.df$Class),]; cell.order <- rownames(annot.df)
# colors
pact.col <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))
class.colors <- ClusterColors(4); names(class.colors) <- sort(unique(model.mvc))
pop.colors <- c('lightgrey', 'black'); names(pop.colors) <- c('AK011', 'AK012')
annot.color <- list('Class' = class.colors, 'Population' = pop.colors)
# plot.dat / breaks
plot.dat <- vip.mat[intersect(unique.ensg, rownames(vip.mat)), cell.order]
plot.dat <- GeneNameConvert(plot.dat, 'human', 'ensg', 'gn')
mat.breaks <- QuantileBreaks(plot.dat, 100)
# PLOT
dev.off()
jpeg('atac/unique-motifs_vip-heatmap.jpg', width = 1200, height = 1200)
pheatmap(plot.dat, main = 'Unique Motif Activity', fontsize = 20,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE,  show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = pact.col(length(mat.breaks) - 1))
dev.off()
###############

## gata activity followup
###############
## modify list
unique.list <- unique.motifs[,2]
unique.list <- unlist(lapply(unique.list, function(x) { unlist(strsplit(x, ',')) }))
unique.list <- unlist(lapply(unique.list, function(x) { unlist(strsplit(x, '::')) }))
unique.list <- gsub("\\(.*", "", unique.list)
## get gatas and convert
gata.list <- unique.list[which(substr(unique.list, 1, 4) == 'GATA')]
gata.ensg <- GeneListConvert(gata.list, 'human', 'gn', 'ensg')
gata.ensg <- intersect(gata.ensg, rownames(vip.mat))
## set classes
ak11.hsc <- intersect(ak11.samps, names(model.mvc[which(model.mvc == 'hsc')]))
ak11.mlp <- intersect(ak11.samps, names(model.mvc[which(model.mvc == 'mlp')]))
ak11.mpp <- intersect(ak11.samps, names(model.mvc[which(model.mvc == 'mpp')]))
ak11.prog <- intersect(ak11.samps, names(model.mvc[which(model.mvc == 'prog')]))
ak12.hsc <- intersect(ak12.samps, names(model.mvc[which(model.mvc == 'hsc')]))
ak12.mlp <- intersect(ak12.samps, names(model.mvc[which(model.mvc == 'mlp')]))
ak12.mpp <- intersect(ak12.samps, names(model.mvc[which(model.mvc == 'mpp')]))
ak12.prog <- intersect(ak12.samps, names(model.mvc[which(model.mvc == 'prog')]))
## viper matrix
gata.vip <- matrix(0L, nrow = length(gata.ensg), ncol = 8)
rownames(gata.vip) <- gata.ensg; colnames(gata.vip) <- c('ak11.hsc', 'ak11.mlp', 'ak11.mpp', 'ak11.prog',
                                                         'ak12.hsc', 'ak12.mlp', 'ak12.mpp', 'ak12.prog')
for (g in gata.ensg) {
  vip.mean.vect <- c(mean(vip.mat[g, ak11.hsc]), mean(vip.mat[g, ak11.mlp]), mean(vip.mat[g, ak11.mpp]), mean(vip.mat[g, ak11.prog]),
                     mean(vip.mat[g, ak12.hsc]), mean(vip.mat[g, ak12.mlp]), mean(vip.mat[g, ak12.mpp]), mean(vip.mat[g, ak12.prog]))
  gata.vip[g,] <- vip.mean.vect
}
## write table
write.table(gata.vip, file = 'atac/unique-gata_class-vip.tsv', sep = '\t', row.names = TRUE, col.names = FALSE, quote = FALSE)
###############



