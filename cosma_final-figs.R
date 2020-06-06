##packages and wd
library(PISCES)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(pheatmap)
library(RColorBrewer)
library(randomForest)
library(plyr)
library(pROC)
library(MASS)
setwd('C://Users/lvlah/linux/ac_lab/cosma/')
source('cosma_functions.R')
set.seed(343)

## set global theme
plot.theme <- theme(
  panel.background = element_rect(fill = 'white', colour = 'white'),
  panel.border = element_rect(color = 'grey', fill = 'transparent'),
  panel.grid = element_line(color = 'grey'),
  plot.title = element_text(hjust = 0.5),
  text = element_text(size = 20)
)
## set global class colors
class.colors <- c(ClusterColors(4), 'grey')
ref.colors <- class.colors; names(ref.colors) <- c('AK001', 'AK002', 'AK003', 'AK004', 'NA')
exp.colors <- class.colors; names(exp.colors) <- c('hsc', 'mpp', 'mlp', 'prog', 'NA')
## ensg names
baz2b.ensg <- 'ENSG00000123636'

### Reference population plots
###############
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vumap.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p01.rds')
marker.list <- c('GATA1', 'GATA2', 'HMGA2', 'BAZ2B', 'BCL11A', 'CDK4', 'HHEX', 'CD53')
## convert to gene symbol names
ref.vip.gs <- Ensemble2GeneName(ref.vip)
## get class vectors
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
train.label <- rep('NA', nrow(ref.umap$layout)); names(train.label) <- rownames(ref.umap$layout)
for (class in unique(train.vect)) {
  class.samps <- names(train.vect)[which(train.vect == class)] 
  train.label[class.samps] <- class
}
## make plot data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'source' = cell.source, 'train' = train.label)
plot.dat <- cbind(plot.dat, t(ref.vip.gs[marker.list,]))
## reference population scatter
jpeg('plots/finalFigs/ref-pops_vip-umap.jpg', width = 1200, height = 1150)
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = source)) + geom_density_2d(color = 'forestgreen') +
  ggtitle('Reference Populations') + labs(color = 'Sample') + plot.theme + scale_color_manual(values = ref.colors)
dev.off()
## reference population selection
jpeg('plots/finalFigs/ref-pops_training-selection_vip-umap.jpg', width = 1200, height = 1150)
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = train)) + scale_color_manual(values = exp.colors) +
  ggtitle('Reference Populations - Training Set Selection') + labs(color = 'Training Set') + plot.theme
dev.off()
## marker plots
for (m in marker.list) {
  jpeg(file = paste('plots/finalFigs/ref-pops_', tolower(m), '-activity_vip-umap.jpg', sep = ''), width = 1200, height = 1150)
  print(ggplot(plot.dat, aes_string(x = 'UMAP1', y = 'UMAP2', color = m)) + geom_point() + geom_density_2d(color = 'forestgreen') +
    ggtitle(m) + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey') + plot.theme)
  dev.off()
}
###############

### AK11/12 plots
###############
exp.cpm <- readRDS('data/exp-samps/exp-pops_cpm.rds')
exp.vip <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
exp.umap <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-1112-vumap.rds'); umap.pref <- 'u'
model.name <- 'dd-p01-pwf'
exp.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-mvc.rds')
exp.percent.mat <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-percent-mat.rds')
## filter for necessary samples
ak11.samps <- colnames(exp.cpm)[which(substr(colnames(exp.cpm), 1, 5) == 'AK011')]
ak12.samps <- colnames(exp.cpm)[which(substr(colnames(exp.cpm), 1, 5) == 'AK012')]
exp.cpm <- exp.cpm[, c(ak11.samps, ak12.samps)]; exp.vip <- exp.vip[, c(ak11.samps, ak12.samps)]
cell.source <- substr(rownames(exp.umap$layout), 1, 5)
## density plots
dense.dat <- data.frame('BAZ2B.a' = exp.vip[baz2b.ensg, c(ak11.samps, ak12.samps)], 
                        'BAZ2B.e' = log(exp.cpm[baz2b.ensg,c(ak11.samps, ak12.samps)], base = 2),
                        'pop' = mapvalues(substr(c(ak11.samps, ak12.samps), 1, 5), c('AK011', 'AK012'), c('Luciferase', 'BAZ2B')))
jpeg('plots/finalFigs/ak1112_baz2b-expression.jpg', width = 1000, height = 1150)
ggplot(dense.dat) + geom_density(aes(x = BAZ2B.e, group = pop, fill = pop), alpha = 0.6) + ggtitle('BAZ2B Expression') +
  xlab('Log2 CPM') + ylab('Density') + plot.theme + theme(legend.title = element_blank())
dev.off()
jpeg('plots/finalFigs/ak1112_baz2b-activity.jpg', width = 1000, height = 1150)
ggplot(dense.dat) + geom_density(aes(x = BAZ2B.a, group = pop, fill = pop), alpha = 0.6) + ggtitle('BAZ2B Activity') +
  xlab('Activity') + ylab('Density') + plot.theme + theme(legend.title = element_blank())
dev.off()
## combined plot
p1 <- ggplot(dense.dat) + geom_density(aes(x = BAZ2B.e, group = pop, fill = pop), alpha = 0.6) + ggtitle('BAZ2B Expression') +
  xlab('Log2 CPM') + ylab('Density') + plot.theme + theme(legend.title = element_blank())
p2 <- ggplot(dense.dat) + geom_density(aes(x = BAZ2B.a, group = pop, fill = pop), alpha = 0.6) + ggtitle('BAZ2B Activity') +
  xlab('Activity') + ylab('Density') + plot.theme + theme(legend.title = element_blank())
dense.plot <- ggarrange(plotlist = list(p1, p2), nrow = 1, ncol = 2)
jpeg('plots/finalFigs/ak1112_baz2b-density.jpg', width = 2000, height = 1150)
print(dense.plot)
dev.off()
## calculate classification entropy
class.entropy.vect <- ClassEntropyVect(exp.percent.mat)
class.alpha <- class.entropy.vect - min(class.entropy.vect); class.alpha <- class.alpha * (0.9 / max(class.alpha))
class.alpha <- 1 - class.alpha
## make plot data
plot.dat <- data.frame('UMAP1' = exp.umap$layout[,1], 'UMAP2' = exp.umap$layout[,2],
                       'source' = cell.source, 'class' = exp.mvc[c(ak11.samps, ak12.samps)],
                       'class.entropy' = class.alpha,
                       'baz2b.p' = exp.vip[baz2b.ensg,], 'baz2b.e' = exp.cpm[baz2b.ensg,])
x.min <- min(exp.umap$layout[,1]) - 0.5; x.max <- max(exp.umap$layout[,1]) + 0.5 
y.min <- min(exp.umap$layout[,2]) - 0.5; y.max <- max(exp.umap$layout[,2]) + 0.5
## pooled classification; w/ and w/out entropy
jpeg(paste('plots/finalFigs/ak1112_', model.name, '_mvc-class_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK011 / AK012 Classification') +
  geom_point(aes(color = class)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class') + plot.theme + scale_color_manual(values = exp.colors)
dev.off()
jpeg(paste('plots/finalFigs/ak1112_', model.name, '_mvc-class-ent_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK011 / AK012 Classification') +
  geom_point(aes(color = class, alpha = class.entropy)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class', alpha = 'Entropy') + 
  plot.theme + scale_color_manual(values = exp.colors)
dev.off()
## AK11 classifications; w/ and w/out entropy
jpeg(paste('plots/finalFigs/ak11_', model.name, '_mvc-class_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak11.samps,], aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK011 Classification') +
  geom_point(aes(color = class)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class') + plot.theme + scale_color_manual(values = exp.colors)
dev.off()
jpeg(paste('plots/finalFigs/ak11_', model.name, '_mvc-class-ent_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak11.samps,], aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK011 Classification') +
  geom_point(aes(color = class, alpha = class.entropy)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class', alpha = 'Entropy') + 
  plot.theme + scale_color_manual(values = exp.colors)
dev.off()
## AK12 classifications; w/ and w/out entropy
jpeg(paste('plots/finalFigs/ak12_', model.name, '_mvc-class_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak12.samps,], aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK012 Classification') +
  geom_point(aes(color = class)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class') + plot.theme + scale_color_manual(values = exp.colors)
dev.off()
jpeg(paste('plots/finalFigs/ak12_', model.name, '_mvc-class-ent_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak12.samps,], aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK012 Classification') +
  geom_point(aes(color = class, alpha = class.entropy)) + geom_density_2d(color = 'forestgreen') + labs(color = 'Class', alpha = 'Entropy') + 
  plot.theme + scale_color_manual(values = exp.colors)
dev.off()
## BAZ2B actiivty scatters
jpeg(paste('plots/finalFigs/ak11_baz2b-pAct_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak11.samps,], aes(UMAP1, UMAP2)) + geom_point(aes(color = baz2b.p)) + xlim(x.min, x.max) + ylim(y.min, y.max) +
  ggtitle('AK011 BAZ2B Activity') + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = mean(plot.dat[ak11.samps, 'baz2b.p'])) + 
  geom_density_2d(color = 'forestgreen') + labs(color = 'BAZ2B') + plot.theme
dev.off()
jpeg(paste('plots/finalFigs/ak12_baz2b-pAct_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak12.samps,], aes(UMAP1, UMAP2)) + geom_point(aes(color = baz2b.p)) + xlim(x.min, x.max) + ylim(y.min, y.max) +
  ggtitle('AK012 BAZ2B Activity') + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = mean(plot.dat[ak12.samps, 'baz2b.p'])) + 
  geom_density_2d(color = 'forestgreen') + labs(color = 'BAZ2B') + plot.theme
dev.off()
## BAZ2B expression scatters
jpeg(paste('plots/finalFigs/ak11_baz2b-gExp_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak11.samps,], aes(UMAP1, UMAP2)) + geom_point(aes(color = baz2b.e)) + xlim(x.min, x.max) + ylim(y.min, y.max) +
  ggtitle('AK011 BAZ2B Expression') + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = mean(plot.dat[ak11.samps, 'baz2b.e'])) + 
  geom_density_2d(color = 'forestgreen') + labs(color = 'BAZ2B') + plot.theme
dev.off()
jpeg(paste('plots/finalFigs/ak12_baz2b-gExp_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
ggplot(plot.dat[ak12.samps,], aes(UMAP1, UMAP2)) + geom_point(aes(color = baz2b.e)) + xlim(x.min, x.max) + ylim(y.min, y.max) +
  ggtitle('AK012 BAZ2B Expression') + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = mean(plot.dat[ak12.samps, 'baz2b.e'])) + 
  geom_density_2d(color = 'forestgreen') + labs(color = 'BAZ2B') + plot.theme
dev.off()
## Vehicle Expression scatter
rfp.gene <- '5LTR_BAZ2B_UbC_rtTA_GFP_3LTR_gene'
luciferase.gene <- 'Luciferase_cDNA_gene'
## new circle plots
CirclePlot(exp.mvc[ak11.samps], exp.percent.mat[ak11.samps,], class.entropy.vect[ak11.samps], plot.title = 'AK011',
           plot.path = paste('plots/finalFigs/ak11_', model.name, '_circle-plot.jpg', sep = ''), exp.colors)
CirclePlot(exp.mvc[ak12.samps], exp.percent.mat[ak12.samps,], class.entropy.vect[ak12.samps], plot.title = 'AK012',
           plot.path = paste('plots/finalFigs/ak12_', model.name, '_circle-plot.jpg', sep = ''), exp.colors)
###############

## ATAC TF Plots
###############
exp.vip <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
exp.umap <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-1112-vumap.rds'); umap.pref <- 'u'
exp.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-mvc.rds')
exp.percent.mat <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-percent-mat.rds')
class.entropy.vect <- ClassEntropyVect(exp.percent.mat)
## filter for necessary samples
ak11.samps <- colnames(exp.vip)[which(substr(colnames(exp.vip), 1, 5) == 'AK011')]
ak12.samps <- colnames(exp.vip)[which(substr(colnames(exp.vip), 1, 5) == 'AK012')]
exp.vip <- exp.vip[, c(ak11.samps, ak12.samps)]
tf.list <- c('FOS', 'GATA2', 'ERG')
## plot.dat
exp.vip.gn <- GeneNameConvert(exp.vip, species = 'human', 'ensg', 'gn')
plot.dat <- data.frame('UMAP1' = exp.umap$layout[,1], 'UMAP2' = exp.umap$layout[,2],
                       'class' = toupper(exp.mvc), 
                       'full.class' = paste(substr(names(exp.mvc), 1, 5), toupper(exp.mvc), sep = '.'))
plot.dat <- cbind(plot.dat, t(exp.vip.gn[tf.list, c(ak11.samps, ak12.samps)]))
x.min <- min(exp.umap$layout[,1]) - 0.5; x.max <- max(exp.umap$layout[,1]) + 0.5 
y.min <- min(exp.umap$layout[,2]) - 0.5; y.max <- max(exp.umap$layout[,2]) + 0.5
## umap plots
for (tf in tf.list) {
  print(tf)
  # ak11 plot
  jpeg(paste('plots/finalFigs/ak11_', tf, '-pAct_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
  print(ggplot(plot.dat[ak11.samps,], aes(UMAP1, UMAP2)) + geom_point(aes_string(color = tf)) + 
    xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle(paste('AK011', tf, 'Activity')) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = 0) + 
    geom_density_2d(color = 'forestgreen') + labs(color = tf) + plot.theme)
  dev.off()
  # ak12 plot
  jpeg(paste('plots/finalFigs/ak12_', tf, '-pAct_vip-', umap.pref, '.jpg', sep = ''), width = 1200, height = 1150)
  print(ggplot(plot.dat[ak12.samps,], aes(UMAP1, UMAP2)) + geom_point(aes_string(color = tf)) + 
    xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle(paste('AK012', tf, 'Activity')) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = 0) + 
    geom_density_2d(color = 'forestgreen') + labs(color = tf) + plot.theme)
  dev.off()
}
## violin plot
for (tf in tf.list) {
  # ak11 plot
  jpeg(paste('plots/finalFigs/ak11-_', tf, '-pAct_vip-violin.jpg', sep = ''), width = 1000, height = 750)
  print(ggplot(plot.dat[ak11.samps,]) + ggtitle(paste('AK011', tf, 'Protein Activity')) + xlab('Classification') + ylab('VIPER NES') + 
          geom_violin(aes_string(x = 'class', y = tf, color = 'class', fill = 'class'), alpha = 0.5) +
          plot.theme + theme(legend.title = element_blank()))
  dev.off()
  # ak12 plot
  jpeg(paste('plots/finalFigs/ak12_', tf, '-pAct_vip-violin.jpg', sep = ''), width = 1000, height = 750)
  print(ggplot(plot.dat[ak12.samps,]) + ggtitle(paste('AK012', tf, 'Protein Activity')) + xlab('Classification') + ylab('VIPER NES') + 
          geom_violin(aes_string(x = 'class', y = tf, color = 'class', fill = 'class'), alpha = 0.5) +
          plot.theme + theme(legend.title = element_blank()))
  dev.off()
  # combined plot w/ separate labels
  jpeg(paste('plots/finalFigs/ak11-12_', tf, '-pAct_vip-violin.jpg', sep = ''), width = 1200, height = 750)
  print(ggplot(plot.dat) + geom_violin(aes_string(x = 'full.class', y = tf, color = 'full.class', fill = 'full.class'), alpha = 0.5) +
    ggtitle(paste(tf, 'Protein Activity')) + xlab('Classification') + ylab('VIPER NES') + 
    plot.theme + theme(legend.title = element_blank()))
  dev.off()
  # combined plot w/ pooled labels
  jpeg(paste('plots/finalFigs/ak11-12_', tf, '-pAct_vip-violin-pooled.jpg', sep = ''), width = 1000, height = 750)
  print(ggplot(plot.dat) + geom_violin(aes_string(x = 'class', y = tf, color = 'class', fill = 'class'), alpha = 0.5) +
          ggtitle(paste(tf, 'Protein Activity')) + xlab('Classification') + ylab('VIPER NES') + 
          plot.theme + theme(legend.title = element_blank()))
  dev.off()
}
## circle plot w/ activity
for (tf in tf.list) {
  # ak11 plot
  CirclePlot(exp.mvc[ak11.samps], exp.percent.mat[ak11.samps,], class.entropy.vect[ak11.samps], 
             plot.title = paste('AK011: ', tf, ' Activity', sep = ''),
             plot.path = paste('plots/finalFigs/ak11_', tf, '-activity_circle-plot.jpg', sep = ''), exp.colors,
             p.vect = plot.dat[ak11.samps, tf], p.name = tf)
  # ak12 plot
  CirclePlot(exp.mvc[ak12.samps], exp.percent.mat[ak12.samps,], class.entropy.vect[ak12.samps], 
             plot.title = paste('AK012: ', tf, ' Activity', sep = ''),
             plot.path = paste('plots/finalFigs/ak12_', tf, '-activity_circle-plot.jpg', sep = ''), exp.colors,
             p.vect = plot.dat[ak12.samps, tf], p.name = tf)
}
###############




