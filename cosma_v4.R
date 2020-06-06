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

### model optimization script
###############
ref.mat <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
ref.labels <- substr(colnames(ref.mat), 1, 5)
ref.labels <- mapvalues(ref.labels, c('AK001', 'AK002', 'AK003', 'AK004'), c('hsc', 'mpp', 'mlp', 'prog'))
## train / test split
print('Performing train / test split...')
tt.split <- TrainTestSplit(ref.mat, ref.labels, train.per = 0.7)
## feature selection
print('Identifying candidate features...')
if (feature.option == 'kw') {
  candidate.features <- kwFeats(tt.split$train.x, tt.split$train.y)
} else if (feature.option == 'pwf') {
  candidate.features <- pwFeats(tt.split$train.x, tt.split$train.y)
} else if (feature.option == 'pwff') {
  candidate.features <- pwFeatsFull(tt.split$train.x, tt.split$train.y)
} else if (feature.option == 'pwffu') {
  candidate.features <- pwFeatsFullUnique(tt.split$train.x, tt.split$train.y)
}  
## featue number optimization
print('Optimizing feature set...')
feat.iter <- length(candidate.features); feat.auroc <- c()
for (fi in 1:feat.iter) {
  print(fi)
  # train model
  model.feats <- candidate.features[[fi]]
  rf.model <- randomForest(x = t(tt.split$train.x[model.feats,,drop=FALSE]), y = tt.split$train.y, ntree = 1000)
  # test model
  mc.auroc <- mcAUROC(rf.model, tt.split$test.x[model.feats,], tt.split$test.y)
  feat.auroc <- c(feat.auroc, mean(mc.auroc))
}
opt.feats <- candidate.features[[which.max(feat.auroc)]]; num.feats <- length(opt.feats)
## mtry optimization
print('Optimizing mtry...')
mtry.iter <- 2:ceiling(sqrt(num.feats))
mtry.auroc <- c()
for (mt in mtry.iter) {
  print(mt)
  # train model
  rf.model <- randomForest(x = t(tt.split$train.x[opt.feats,]), y = tt.split$train.y, ntree = 1000, mtry = mt)
  # test model
  mc.auroc <- mcAUROC(rf.model, tt.split$test.x[opt.feats,], tt.split$test.y)
  mtry.auroc <- c(mtry.auroc, mean(mc.auroc))
}
opt.mtry <- mtry.iter[which.max(mtry.auroc)]
## train final model
print('Training final model...')
final.model <- randomForest(x = t(tt.split$train.x[opt.feats,]), y = tt.split$train.y, ntree = 10000, mtry = opt.mtry)
## Youden's method for ROC
print('Optimizing class thresholds...')
test.predict <- predict(final.model, t(tt.split$test.x[opt.feats,]), predict.all = TRUE)
class.labels <- levels(tt.split$test.y); ntree <- ncol(test.predict$individual)
yt.thresh <- list()
for (cl in class.labels) {
  print(cl)
  # get class ROC
  class.vect <- as.character(tt.split$test.y); class.vect[which(class.vect != cl)] <- 'other'
  class.percent <- apply(test.predict$individual, 1, function(x) {
    return(length(which(x == cl)) / ntree)
  })
  class.roc <- roc(class.vect, class.percent)
  # get youden's threshold; create vector of youden's index, find max as threshold
  class.yi <- class.roc$sensitivities + class.roc$specificities - 1; max.yi <- which.max(class.yi)
  class.yt <- class.roc$thresholds[max.yi]
  yt.thresh[[cl]] <- class.yt
}
## save final model
print('Saving final model...')
model.save.obj <- list('model' = final.model, 'model.feats' = opt.feats, 'model.yt' = yt.thresh,
                      'feat.auroc' = feat.auroc, 'mtry.auroc' = mtry.auroc, 'train.samps' = colnames(tt.split$train.x))
saveRDS(model.save.obj, file = 'model-results/k10rnv_kw_model-obj.rds')
###############

### model testing script
###############
exp.vip <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
model.obj <- readRDS('model-results/model-objs/k10rnv_leR-v2_pwffu.rds')
exp.labels <- substr(colnames(exp.vip), 1, 5)
f.name <- 'k10rnv-leR-pwffu-2_lnorm'
## generate predictiosn and percentage matrix
exp.predict <- predict(model.obj$model, t(exp.vip[model.obj$model.feats,]), predict.all = TRUE)
percent.mat <- MakePercentMat(exp.predict)
class.labels <- levels(exp.predict$aggregate)
conf.table <- table(exp.predict$aggregate, exp.labels)
conf.table <- round(t(t(conf.table) / colSums(conf.table)), digits = 3)
write.csv(conf.table, file = paste('model-results/conf-tables/', f.name, '_ctable.csv', sep = ''))
## yi classify
yi.mat <- percent.mat
for (cl in class.labels) {
  yi.vect <- rep(0, nrow(yi.mat))
  yi.vect[which(yi.mat[,cl] >= model.obj$model.yt[cl])] <- 1
  yi.mat[,cl] <- yi.vect
}
yi.class <- apply(yi.mat, 1, function(x) {
  if (sum(x) != 1) {
    return('unknown')
  } else {
    return(colnames(yi.mat)[which.max(x)])
  }
})
yi.conf.table <- table(yi.class, as.factor(exp.labels))
yi.conf.table <- round(t(t(yi.conf.table) / colSums(yi.conf.table)), digits = 3)
write.csv(yi.conf.table, file = paste('model-results/conf-tables/', f.name, '_yi-ctable.csv', sep = ''))
## get fisher
yi.hsc.0506 <- yi.class[which(substr(names(yi.class), 1, 5) %in% c('AK005', 'AK006'))]
yi.hsc.0506[which(yi.hsc.0506 != 'hsc')] <- 'other'
fisher.test(as.factor(yi.hsc.0506), as.factor(substr(names(yi.hsc.0506), 1, 5)))
## get fisher for 11/12
yi.hsc.1112 <- yi.class[which(substr(names(yi.class), 1, 5) %in% c('AK011', 'AK012'))]
yi.hsc.1112[which(yi.hsc.1112 != 'hsc')] <- 'other'
fisher.test(as.factor(yi.hsc.1112), as.factor(substr(names(yi.hsc.1112), 1, 5)))
###############

### comprehensive model testing
###############
## load e56 matrices and subset for memory space
e56.rnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-rnorm-vip.rds')
e56.rnorm <- e56.rnorm[, which(substr(colnames(e56.rnorm), 1, 5) %in% c('AK005', 'AK006'))]
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
e56.lnorm <- e56.lnorm[, which(substr(colnames(e56.lnorm), 1, 5) %in% c('AK005', 'AK006'))]
e56.rnorm.labels <- substr(colnames(e56.rnorm), 1, 5)
e56.lnorm.labels <- substr(colnames(e56.lnorm), 1, 5)
## load network-specific rnorm / lnorm
exp.rnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-rnorm-vip.rds')
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
rnorm.labels <- substr(colnames(exp.rnorm), 1, 5)
lnorm.labels <- substr(colnames(exp.lnorm), 1, 5)
## load model and set name
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p0075-pwffu.rds')
model.name <- 'k10rnv_dd-p0075_pwffu'
name.pref <- paste('model-results/conf-tables/', model.name, '/', model.name, sep = '')
## rnorm test
exp.pred <- predict(model.obj$model, t(exp.rnorm[model.obj$model.feats,]), predict.all = TRUE)
percent.mat <- MakePercentMat(exp.pred)
yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
WriteConfTable(exp.pred$aggregate, rnorm.labels, out.file = paste(name.pref, '_rnorm-mvc.csv', sep = ''))
WriteConfTable(yi.pred, rnorm.labels, out.file = paste(name.pref, '_rnorm-yi.csv', sep = ''))
## lnorm test
exp.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
percent.mat <- MakePercentMat(exp.pred)
yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
WriteConfTable(exp.pred$aggregate, lnorm.labels, out.file = paste(name.pref, '_lnorm-mvc.csv', sep = ''))
WriteConfTable(yi.pred, lnorm.labels, out.file = paste(name.pref, '_lnorm-yi.csv', sep = ''))
## e56 rnorm test
exp.pred <- predict(model.obj$model, t(e56.rnorm[model.obj$model.feats,]), predict.all = TRUE)
percent.mat <- MakePercentMat(exp.pred)
yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
WriteConfTable(exp.pred$aggregate, e56.rnorm.labels, out.file = paste(name.pref, '_e56-rnorm-mvc.csv', sep = ''))
WriteConfTable(yi.pred, e56.rnorm.labels, out.file = paste(name.pref, '_e56-rnorm-yi.csv', sep = ''))
## e56 lnorm test
exp.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
percent.mat <- MakePercentMat(exp.pred)
yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
WriteConfTable(exp.pred$aggregate, e56.lnorm.labels, out.file = paste(name.pref, '_e56-lnorm-mvc.csv', sep = ''))
WriteConfTable(yi.pred, e56.lnorm.labels, out.file = paste(name.pref, '_e56-lnorm-yi.csv', sep = ''))
###############

### generate classification objects and save
###############
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p01-pwf.rds'); model.name <- 'k10-e56-dd-p01-pwf'
## load 56 and classify
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
exp5.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK006')]
e56.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
e56.mvc <- e56.pred$aggregate[c(exp5.samps, exp6.samps)]
e56.percent.mat <- MakePercentMat(e56.pred)
## load 1112 and classify
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
exp11.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK012')]
exp.lnorm <- exp.lnorm[,c(exp11.samps, exp12.samps)]
k10rnv.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
e1112.mvc <- k10rnv.pred$aggregate[c(exp11.samps, exp12.samps)]
e1112.percent.mat <- MakePercentMat(k10rnv.pred); e1112.percent.mat <- e1112.percent.mat[c(exp11.samps, exp12.samps),]
## save classification objects
saveRDS(e56.mvc, file = paste('model-results/classification-objects/', model.name, '_e56-mvc.rds', sep = ''))
saveRDS(e56.percent.mat, file = paste('model-results/classification-objects/', model.name, '_e56-percent-mat.rds', sep = ''))
saveRDS(e1112.mvc, file = paste('model-results/classification-objects/', model.name, '_e1112-mvc.rds', sep = ''))
saveRDS(e1112.percent.mat, file = paste('model-results/classification-objects/', model.name, '_e1112-percent-mat.rds', sep = ''))
###############

### model based UMAP
###############
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p02-pwf.rds'); model.name <- 'k10-e56-dd-p02-pwf'
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
exp5.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK006')]
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
exp11.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK012')]
## e56 umap
e56.feature.umap <- CustomUMAP(e56.lnorm[model.obj$model.feats, c(exp5.samps, exp6.samps)])
saveRDS(e56.feature.umap, file = paste('model-results/classification-objects/', model.name, '_e56-feat-umap.rds', sep = ''))
## e1112 umap
e1112.feature.umap <- CustomUMAP(exp.lnorm[model.obj$model.feats, c(exp11.samps, exp12.samps)])
saveRDS(e1112.feature.umap, file = paste('model-results/classification-objects/', model.name, '_e1112-feat-umap.rds', sep = ''))
###############

### mwkmeans clustering
###############
vip.mat <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
vip.dist <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vdist.rds')
clust.vect <- readRDS('clusterings/ref-pops_k10rnv-louvain.rds')
opt.clust <- clust.vect$clusterings[[which.max(clust.vect$sils)]]
## get cluster centers
oc.table <- summary(opt.clust)
oc.centers <- list()
for (k in names(oc.table)) {
  clust.samples <- names(opt.clust)[which(opt.clust == k)]
  oc.centers[[k]] <- rowMeans(vip.mat[, clust.samples])
}
center.mat <- do.call(cbind, oc.centers)
## generate clustering
mwk.clust <- Cluster(vip.mat, center.mat)
saveRDS(mwk.clust, file = 'clusterings/ref-pops_kl0rnv-louvain-mwk.rds')
## plot

###############

### population specific scatters
###############
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
cell.source <- substr(colnames(ref.vip), 1, 5)
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vumap.rds')
## make data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source)
## make plots
plot.title <- 'K10RNV'
plot.list <- list()
for (cl in unique(plot.dat$population)) {
  s.plot.dat <- plot.dat[which(plot.dat$population == cl),]
  plot.list[[cl]] <- ggplot(s.plot.dat, aes(UMAP1, UMAP2)) + geom_point(color = 'darkgrey') + 
    geom_density_2d(color = 'forestgreen') + ggtitle(cl)
}
source.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 2)
jpeg(file = 'plots/population-scatters/k10rnv_pops.jpg', width = 1000, height = 1100)
print(annotate_figure(source.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### plot reference VIPER results
###############
k10rnv.mat <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
bcnv.mat <- readRDS('data/viper-mats/bcnv/ref-pops_bcnv.rds')
baz2b.ensg <- 'ENSG00000123636'
## make data frames
p.dat.1 <- data.frame('BAZ2B' = k10rnv.mat[baz2b.ensg,], 'population' = substr(colnames(k10rnv.mat), 1, 5))
p.dat.2 <- data.frame('BAZ2B' = bcnv.mat[baz2b.ensg,], 'population' = substr(colnames(bcnv.mat), 1, 5))
## plot
plot.title <- 'Referene Population BAZ2B Activity'
p1 <- ggplot(p.dat.1) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('K10RNV') + xlab('BAZ2B Activity') + ylab('Density')
p2 <- ggplot(p.dat.2) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('BCNV') + xlab('BAZ2B Activity') + ylab('Density')
clust.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = 'plots/baz2b-density/ref-pops_baz2b.jpg', width = 1500, height = 750)
print(annotate_figure(clust.plot, top = text_grob(plot.title, family = 'TT Arial', size = 16)))
dev.off()
###############

### plot exp VIPER results
###############
k10rnv.mat <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-rnorm-vip.rds')
bcnv.mat <- readRDS('data/viper-mats/bcnv/exp-pops_bcn-rnorm-vip.rds')
baz2b.ensg <- 'ENSG00000123636'
## make data frames
p.dat.1 <- data.frame('BAZ2B' = k10rnv.mat[baz2b.ensg,], 'population' = substr(colnames(k10rnv.mat), 1, 5))
p.dat.2 <- data.frame('BAZ2B' = bcnv.mat[baz2b.ensg,], 'population' = substr(colnames(bcnv.mat), 1, 5))
## plot
plot.title <- 'Experiment BAZ2B Activity: RNORM'
p1 <- ggplot(p.dat.1) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('K10RNV') + xlab('BAZ2B Activity') + ylab('Density')
p2 <- ggplot(p.dat.2) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('BCNV') + xlab('BAZ2B Activity') + ylab('Density')
clust.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = 'plots/baz2b-density/exp-pops-rnorm_baz2b.jpg', width = 1500, height = 750)
print(annotate_figure(clust.plot, top = text_grob(plot.title, family = 'TT Arial', size = 16)))
dev.off()
###############

### e56 net results
###############
rnorm.mat <- readRDS('data/viper-mats/e56nv/e56_ep56n-rnorm-vip.rds')
rnorm.mat <- rnorm.mat[,which(substr(colnames(rnorm.mat), 1, 5) %in% c('AK005', 'AK006'))]
lnorm.mat <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
lnorm.mat <- lnorm.mat[,which(substr(colnames(lnorm.mat), 1, 5) %in% c('AK005', 'AK006'))]
baz2b.ensg <- 'ENSG00000123636'
## make data frames
p.dat.1 <- data.frame('BAZ2B' = rnorm.mat[baz2b.ensg,], 'population' = substr(colnames(rnorm.mat), 1, 5))
p.dat.2 <- data.frame('BAZ2B' = lnorm.mat[baz2b.ensg,], 'population' = substr(colnames(lnorm.mat), 1, 5))
## plot
plot.title <- 'AK005/6 Network: BAZ2B'
p1 <- ggplot(p.dat.1) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('RNORM') + xlab('BAZ2B Activity') + ylab('Density')
p2 <- ggplot(p.dat.2) + geom_density(aes(x = BAZ2B, group = population, fill = population), alpha = 0.6) +
  ggtitle('LNORM') + xlab('BAZ2B Activity') + ylab('Density')
clust.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = 'plots/baz2b-density/e56nv_baz2b.jpg', width = 1500, height = 750)
print(annotate_figure(clust.plot, top = text_grob(plot.title, family = 'TT Arial', size = 16)))
dev.off()
###############

### entropy plots
###############
ref.umap <- readRDS('data/ref-samps/ref-pops_cpm-pca-umap.rds')
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
ref.ent <- readRDS('entropy_analysis/ref-pops_cpm-sce1.rds')
## make data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source, 'entropy' = ref.ent)
## make plots
plot.title <- 'CPM SCE1'
p1 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = population)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = entropy)) + geom_density_2d(color = 'forestgreen') +
  scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat$entropy))
p3 <- ggplot(plot.dat) + geom_density(aes(x = entropy, group = population, fill = population), alpha = 0.6) +
  xlab('Single Cell Entropy') + ylab('Density')
ent.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = 'plots/entropy/cpm_sce1.jpg', width = 1800, height = 750)
print(annotate_figure(ent.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### cluster plots
###############
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vpca-umap.rds')
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
ref.clust <- readRDS('clusterings/ref-pops_k10rnv-pca30-louvain.rds')
opt.clust <- ref.clust$clusterings[[which.max(ref.clust$sils)]]
## make data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source, 'cluster' = opt.clust)
## plotting
plot.title <- 'k10rnv: PCA-Louvain Clustering'
p1 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = population)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = cluster)) + geom_density_2d(color = 'forestgreen') 
clust.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = 'plots/clusterings/ref-pops_bcnv-pca-louvain.jpg', width = 1500, height = 750)
print(annotate_figure(clust.plot, top = text_grob(plot.title, family = 'TT Arial', size = 16)))
dev.off()
###############

### cluster-based refinement
###############
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vpca-umap.rds')
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
ref.clust <- readRDS('clusterings/ref-pops_k10rnv-pca30-louvain.rds')
opt.clust <- ref.clust$clusterings[[which.max(ref.clust$sils)]]
ref.ent <- readRDS('entropy_analysis/ref-pops_k10rnv-sce2.rds')
## plot data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source, 'cluster' = opt.clust, 'entropy' = ref.ent)
## make plots
plot.title <- 'K10RNV: PCA + Louvain + SCE2 Reference Refinement'
p1 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = population)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = entropy)) + geom_density_2d(color = 'forestgreen') +
  scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat$entropy))
p3 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = cluster)) + geom_density_2d(color = 'forestgreen') 
ent.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = 'plots/model-refinement/k10rnv-pca-louvain-sce2.jpg', width = 1800, height = 750)
print(annotate_figure(ent.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
## aggregate cluster entropies
cluster.df <- data.frame('cluster' = opt.clust, 'entropy' = ref.ent)
aggregate(cluster.df$entropy, list(cluster.df$cluster), mean)
## make training labels
clust.map <- list('hsc' = 6, 'mpp' = 7, 'mlp' = 2, 'prog' = 10)
sample.names <- c(); sample.labels <- c()
for (cm in names(clust.map)) {
  # get samples
  cm.num <- clust.map[[cm]]
  cm.samps <- which(opt.clust == cm.num)
  # add to lists
  sample.names <- c(sample.names, names(opt.clust)[cm.samps])
  sample.labels <- c(sample.labels, rep(cm, length(cm.samps)))
}
names(sample.labels) <- sample.names
## sanity check and save
table(sample.labels, substr(names(sample.labels), 1, 5))
saveRDS(sample.labels, file = 'model-results/refined-labels/k10rnv_pca-louvain-ent-optimized_v2.rds')
###############

### cluster-based refinement w/ silhouette
###############
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vpca-umap.rds')
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
ref.clust <- readRDS('clusterings/ref-pops_k10rnv-pca30-louvain.rds')
opt.clust <- ref.clust$clusterings[[which.max(ref.clust$sils)]]
ref.ent <- readRDS('entropy_analysis/ref-pops_k10rnv-sce2.rds')
## calculate silhouette score
ref.dist <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vpca-dist.rds')
clust.sil <- silhouette(as.numeric(opt.clust), ref.dist)
sil.vect <- clust.sil[,3]; names(sil.vect) <- names(opt.clust)
## plot data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source, 'cluster' = opt.clust, 'entropy' = ref.ent, 'silhouette' = sil.vect)
## make plots
plot.title <- 'K10RNV: PCA + Louvain + Silhouette + SCE2 Reference Refinement'
p1 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = population)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = entropy)) + geom_density_2d(color = 'forestgreen') +
  scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat$entropy))
p3 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = cluster, alpha = silhouette)) + geom_density_2d(color = 'forestgreen') +
  scale_alpha('silhouette', range = c(0.1,1))
ent.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = 'plots/model-refinement/k10rnv-pca=louvain-sil-sce2.jpg', width = 1800, height = 750)
print(annotate_figure(ent.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
## aggregate cluster entropies
cluster.df <- data.frame('cluster' = opt.clust, 'entropy' = ref.ent)
aggregate(cluster.df$entropy, list(cluster.df$cluster), mean)
## make training labels
clust.map <- list('hsc' = 6, 'mpp' = 7, 'mlp' = 2, 'prog' = 10)
sample.names <- c(); sample.labels <- c()
for (cm in names(clust.map)) {
  # get samples
  cm.num <- clust.map[[cm]]
  cm.samps <- names(opt.clust)[which(opt.clust == cm.num)]
  cm.sil <- sort(sil.vect[cm.samps], decreasing = TRUE)
  cm.samps <- names(cm.sil)[1:300]
  # add to lists
  sample.names <- c(sample.names, cm.samps)
  sample.labels <- c(sample.labels, rep(cm, length(cm.samps)))
}
names(sample.labels) <- sample.names
## sanity check and save
table(sample.labels, substr(names(sample.labels), 1, 5))
saveRDS(sample.labels, file = 'model-results/refined-labels/k10rnv_pca-louvain-silhouette-ent-optimized.rds')
###############

### lymphoid-myeloid plots
###############
ref.mat <- readRDS('data/viper-mats/bcnv/ref-pops_bcnv.rds')
ref.umap <- readRDS('data/viper-mats/bcnv/ref-pops_bcnv-umap.rds')
ref.mat.gn <- Ensemble2GeneName(ref.mat)
myeloid.markers <- c('CEBPA, GATA2, MYBL2, CEBPB')
lymphoid.markers <- c('HHEX', 'BCL11A', 'BCL6', 'CREB1')
lm.markers <- intersect(rownames(ref.mat.gn), c(myeloid.markers, lymphoid.markers))
## plot data frame
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2])
## make plot
plot.title <- 'BCNV: ML Markers'
plot.list <- list()
for (lm in lymphoid.markers) {
  plot.dat[[lm]] <- ref.mat.gn[lm,]
  plot.list[[lm]] <- ggplot(plot.dat, aes_string(x = 'UMAP1', y = 'UMAP2', color = lm)) + geom_point() + ggtitle(lm) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat[,lm]))
}
ml.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 2)
jpeg(file = 'plots/model-refinement/bcnv-mlMarkers.jpg', width = 1400, height = 1500)
print(annotate_figure(ml.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### hsc marker plots
###############
ref.mat <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vumap.rds')
ref.mat.gn <- Ensemble2GeneName(ref.mat)
ref.ent <- readRDS('entropy_analysis/ref-pops_k10rnv-sce2.rds')
cell.source <- substr(colnames(ref.mat.gn), 1, 5)
## set markers
hsc.markers <- c('CDC42', 'SETD2', 'TCF7', 'ARID4B', 'CHD1', 'KLF7', 'MAPKAPK5', 'ZFP36LI', 'HMGA2',
                 'SNAI2', 'DEK', 'SMARCA2', 'E4F1', 'SCMH1', 'PDK1', 'PARP1', 'RHEB', 'AH1L', 'GATA2',
                 'GAB2', 'ARID4A')
hsc.markers <- sort(intersect(hsc.markers, rownames(ref.mat.gn)))
## plot data frame, entropy, and population
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'population' = cell.source, 'entropy' = ref.ent)
plot.list <- list()
plot.list[['population']] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = population)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('Population')
plot.list[['entropy']] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = entropy)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('Entropy') + 
  scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat$entropy))
## marker plots
for (i in hsc.markers) {
  plot.dat[[i]] <- ref.mat.gn[i,]
  plot.list[[i]] <- ggplot(plot.dat, aes_string(x = 'UMAP1', y = 'UMAP2', color = i)) + geom_point() + ggtitle(i) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat[,i]))
}
## generate final plot
plot.title <- 'K10RNV: HSC Markers'
hsc.plot <- ggarrange(plotlist = plot.list, nrow = 4, ncol = 5)
jpeg(file = 'plots/model-refinement/k10rnv_umap-hsc-markers.jpg', width = 2000, height = 1800)
print(annotate_figure(hsc.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### differntial density
###############
ref.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vumap.rds')
cell.source <- substr(rownames(ref.umap$layout), 1, 5)
## overall density
num.grids <- 1000
x.lims <- range(ref.umap$layout[,1]); y.lims <- range(ref.umap$layout[,2])
c.lims <- c(x.lims, y.lims)
total.density <- kde2d(ref.umap$layout[,1], ref.umap$layout[,2], n = num.grids, lims = c.lims)
total.dm <- total.density$z / sum(total.density$z)
x.ss <- total.density$x[2] - total.density$x[1]
y.ss <- total.density$y[2] - total.density$y[1]
## get population samples
ak1.samps <- which(cell.source == 'AK001')
ak2.samps <- which(cell.source == 'AK002')
ak3.samps <- which(cell.source == 'AK003')
ak4.samps <- which(cell.source == 'AK004')
## get population densities
ak1.density <- kde2d(ref.umap$layout[ak1.samps, 1], ref.umap$layout[ak1.samps, 2], n = num.grids, lims = c.lims)
ak2.density <- kde2d(ref.umap$layout[ak2.samps, 1], ref.umap$layout[ak2.samps, 2], n = num.grids, lims = c.lims)
ak3.density <- kde2d(ref.umap$layout[ak3.samps, 1], ref.umap$layout[ak3.samps, 2], n = num.grids, lims = c.lims)
ak4.density <- kde2d(ref.umap$layout[ak4.samps, 1], ref.umap$layout[ak4.samps, 2], n = num.grids, lims = c.lims)
## save density objects
dense.obj <- list('total' = total.density, 'ak1' = ak1.density, 'ak2' = ak2.density, 'ak3' = ak3.density, 'ak4' = ak4.density)
saveRDS(dense.obj, file = 'model-results/k10rnv-umap-density_n1000.rds')
## normalize densities
ak1.dm <- (ak1.density$z / sum(ak1.density$z)) - total.dm
ak2.dm <- (ak2.density$z / sum(ak2.density$z)) - total.dm
ak3.dm <- (ak3.density$z / sum(ak3.density$z)) - total.dm
ak4.dm <- (ak4.density$z / sum(ak4.density$z)) - total.dm
## plot densities of the density delta
plot.dat <- data.frame('dif.density' = c(ak1.dm, ak2.dm, ak3.dm, ak4.dm),
                       'Sample' = c(rep('AK001', num.grids**2), rep('AK002', num.grids**2), rep('AK003', num.grids**2), rep('AK004', num.grids**2)))
ggplot(plot.dat) + geom_density(aes(x = dif.density, group = Sample, fill = Sample), alpha = 0.6) + 
  xlab('UMAP Density Delta') + ylab('Density') + ggtitle('K10RNV Differential Density')
## gets all cells in grids with density delta greater than the threshold
DThreshCells <- function(layout, norm.density, full.density, dt.per = 0.01) {
  # calcualte density thresh
  dense.thresh <- sort(norm.density, decreasing = TRUE)[length(norm.density) * dt.per]
  # select coordinates
  max.coords <- which(norm.density >= dense.thresh, arr.ind = T)
  dense.samps <- c()
  for (i in 1:nrow(max.coords)) {
    # get coordinates
    grid.coords <- max.coords[i,]
    x.coord <- full.density$x[grid.coords[1]]
    y.coord <- full.density$y[grid.coords[2]]
    # get cells
    cell.samps <- GetCells(ref.umap$layout, xmin = (x.coord - x.ss*0.5), xmax = (x.coord + x.ss*0.5),
                           ymin = (y.coord - y.ss*0.5), ymax = (y.coord + y.ss*0.5))
    # add to list
    dense.samps <- c(dense.samps, cell.samps)
  }
  return(unique(dense.samps))
}
## cell fetch function
GetCells <- function(layout, xmin, xmax, ymin, ymax) {
  x.match <- intersect(which(layout[,1] >= xmin), which(layout[,1] <= xmax))
  y.match <- intersect(which(layout[,2] >= ymin), which(layout[,2] <= ymax))
  ret.cells <- rownames(layout)[intersect(x.match, y.match)]
  return(ret.cells)
}
## funciton to get cells from maixmum density
DensityCells <- function(layout, norm.density, full.density, num.grids = 50) {
  max.coords <- which(norm.density >= sort(norm.density, decreasing = T)[num.grids], arr.ind = T)
  dense.samps <- c()
  for (i in 1:nrow(max.coords)) {
    # get coordinates
    grid.coords <- max.coords[i,]
    x.coord <- full.density$x[grid.coords[1]]
    y.coord <- full.density$y[grid.coords[2]]
    # get cells
    cell.samps <- GetCells(ref.umap$layout, xmin = (x.coord - x.ss*0.5), xmax = (x.coord + x.ss*0.5),
                           ymin = (y.coord - y.ss*0.5), ymax = (y.coord + y.ss*0.5))
    # add to list
    dense.samps <- c(dense.samps, cell.samps)
  }
  return(unique(dense.samps))
}
## select cells with density specific percentage
dt.per <- 0.0075
ak1.cells <- DThreshCells(ref.umap$layout, ak1.dm, ak1.density, dt.per)
ak2.cells <- DThreshCells(ref.umap$layout, ak2.dm, ak2.density, dt.per)
ak3.cells <- DThreshCells(ref.umap$layout, ak3.dm, ak3.density, dt.per)
ak4.cells <- DThreshCells(ref.umap$layout, ak4.dm, ak4.density, dt.per)
sample.labels <- c(rep('hsc', length(ak1.cells)), rep('mpp', length(ak2.cells)),
                   rep('mlp', length(ak3.cells)), rep('prog', length(ak4.cells)))
names(sample.labels) <- c(ak1.cells, ak2.cells, ak3.cells, ak4.cells)
saveRDS(sample.labels, file = 'model-results/refined-labels/k10rnv_dif-density-labels_p0075.rds')
## select cells purely based on the number of grids 
#num.grids <- 25
#ak1.cells <- DensityCells(ref.umap$layout, ak1.dm, ak1.density, num.grids = num.grids)
#ak2.cells <- DensityCells(ref.umap$layout, ak2.dm, ak2.density, num.grids = num.grids)
#ak3.cells <- DensityCells(ref.umap$layout, ak3.dm, ak3.density, num.grids = num.grids)
#ak4.cells <- DensityCells(ref.umap$layout, ak4.dm, ak4.density, num.grids = num.grids)
#sample.labels <- c(rep('hsc', length(ak1.cells)), rep('mpp', length(ak2.cells)),
#                   rep('mlp', length(ak3.cells)), rep('prog', length(ak4.cells)))
#names(sample.labels) <- c(ak1.cells, ak2.cells, ak3.cells, ak4.cells)
## plots
cell.labels <- rep('NA', nrow(ref.umap$layout)); names(cell.labels) <- rownames(ref.umap$layout)
for (u in unique(sample.labels)) {
  cell.labels[names(sample.labels)[which(sample.labels == u)]] <- u
}
plot.dat <- data.frame('UMAP1' = ref.umap$layout[,1], 'UMAP2' = ref.umap$layout[,2],
                       'label' = cell.labels)
class.cols <- c(ClusterColors(4), 'grey'); names(class.cols) <- c('hsc', 'mpp', 'mlp', 'prog', 'NA')
jpeg(file = 'plots/model-refinement/k10rnv-umap-dif-density-p0075.jpg', width = 1400, height = 1500)
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = label)) + ggtitle('Differential Density Refinement') +
  scale_color_manual(values = class.cols)
dev.off()
###############

### e1112 result plotting
###############
e1112.umap <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e1112-feat-umap.rds'); umap.pref <- 'fu'
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
exp11.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK012')]
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p02-pwf.rds'); model.name <- 'k10-dd-p02-pwf'
## classify
k10rnv.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
e1112.mvc <- k10rnv.pred$aggregate[c(exp11.samps, exp12.samps)]
e1112.percent.mat <- MakePercentMat(k10rnv.pred); e1112.percent.mat <- e1112.percent.mat[c(exp11.samps, exp12.samps),]
class.colors <- ClusterColors(4); names(class.colors) <- colnames(e1112.percent.mat)
## generate color vectors
class.thresh.col.vect <- ClassThreshColVect(e1112.percent.mat, class.colors, class.thresh = 0.2)
class.entropy.vect <- ClassEntropyVect(e1112.percent.mat)
class.alpha <- class.entropy.vect - min(class.entropy.vect); class.alpha <- class.alpha * (0.9 / max(class.alpha))
class.alpha <- 1 - class.alpha
## generate data frame
ak11.class <- as.vector(e1112.mvc); names(ak11.class) <- names(e1112.mvc); ak11.class[exp12.samps] <- 'AK012'
ak12.class <- as.vector(e1112.mvc); names(ak12.class) <- names(e1112.mvc); ak12.class[exp11.samps] <- 'AK011'
full.class.colors <- class.colors; full.class.colors[['AK011']] <- 'grey'; full.class.colors[['AK012']] <- 'grey'
plot.dat <- data.frame('UMAP1' = e1112.umap$layout[,1], 'UMAP2' = e1112.umap$layout[,2],
                       'class' = e1112.mvc, 'colors' = class.thresh.col.vect, 'alpha' = I(class.alpha),
                       'entropy' = class.entropy.vect, 'sample' = substr(rownames(e1112.umap$layout), 1, 5),
                       'ak11class' = ak11.class, 'ak12class' = ak12.class)
## mvc class plot
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = ak11class, alpha = sample)) +
  scale_alpha_discrete(range = c(1, 0.1)) + scale_color_manual(values = full.class.colors) +
  ggtitle('AK011') + geom_density_2d(color = 'forestgreen')
p2 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = ak12class, alpha = sample)) +
  scale_alpha_discrete(range = c(0.1, 1)) + scale_color_manual(values = full.class.colors) + 
  ggtitle('AK012') + geom_density_2d(color = 'forestgreen')
e1112.mvc.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
print(annotate_figure(e1112.mvc.plot, top = text_grob('AK011 / AK012 Classification', size = 16)))
## sample specific data frames
ak11.alpha <- class.alpha[exp11.samps]; ak12.alpha <- class.alpha[exp12.samps]
ak11.dat <- data.frame('UMAP1' = e1112.umap$layout[exp11.samps, 1], 'UMAP2' = e1112.umap$layout[exp11.samps, 2],
                       'class' = e1112.mvc[exp11.samps], 'alpha' = I(ak11.alpha),
                       'class.mix' = class.thresh.col.vect[exp11.samps])
ak12.dat <- data.frame('UMAP1' = e1112.umap$layout[exp12.samps, 1], 'UMAP2' = e1112.umap$layout[exp12.samps, 2],
                       'class' = e1112.mvc[exp12.samps], 'alpha' = I(ak12.alpha),
                       'class.mix' = class.thresh.col.vect[exp12.samps])
ak11.dat <- cbind(ak11.dat, e1112.percent.mat[exp11.samps,])
ak12.dat <- cbind(ak12.dat, e1112.percent.mat[exp12.samps,])
x.min <- min(e1112.umap$layout[,1]) - 0.5; x.max <- max(e1112.umap$layout[,1]) + 0.5 
y.min <- min(e1112.umap$layout[,2]) - 0.5; y.max <- max(e1112.umap$layout[,2]) + 0.5
## separate mvc plots
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('Combined') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
p2 <- ggplot(ak11.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('AK11') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
p3 <- ggplot(ak12.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('AK12') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
e1112.mvc.sep <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e1112-mvc-sep_', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e1112.mvc.sep, top = text_grob('AK011 / AK012 Classification', size = 16)))
dev.off()
## separate mvc plots w/ entropy alpha
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('Pooled') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(ak11.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK11') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen') 
p3 <- ggplot(ak12.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK12') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen')
e1112.mvc.sep.ent <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e1112-mvc-sep-ent_', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e1112.mvc.sep.ent, top = text_grob('AK011 / AK012 Classification', size = 16)))
dev.off()
## separate mvc plots w/ color mixing
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('Pooled') +
  geom_point(aes(color = colors)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
p2 <- ggplot(ak11.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK11') +
  geom_point(aes(color = class.mix)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
p3 <- ggplot(ak12.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK12') +
  geom_point(aes(color = class.mix)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
e1112.mvc.sep.col.mix <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e1112-mvc-sep-col-mix_', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e1112.mvc.sep.col.mix, top = text_grob('AK011 / AK012 Classification', size = 16)))
dev.off()
## eight pannel plot
plot.list <- list()
for (i in names(class.colors)) {
  # AK11 plot
  plot.obj <- ggplot(ak11.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(ak11.dat[[i]])) + 
    ggtitle(paste('AK11: ', toupper(i), ' Classification', sep = '')) + xlim(x.min, x.max) + ylim(y.min, y.max)
  plot.list[[paste('AK11', i, sep = '.')]] <- plot.obj
  # AK12 plot
  plot.obj <- ggplot(ak12.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(ak12.dat[[i]])) + 
    ggtitle(paste('AK12: ', toupper(i), ' Classification', sep = '')) + xlim(x.min, x.max) + ylim(y.min, y.max)
  plot.list[[paste('AK12', i, sep = '.')]] <- plot.obj
}
plot.list <- plot.list[sort(names(plot.list))]
e1112.eight.panel <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 4)
jpeg(file = paste('plots/model-results/e1112-eight-panel', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 2000, height = 1100)
print(annotate_figure(e1112.eight.panel, top = text_grob('AK011 / AK012 Classification', family = 'Arial', size = 16)))
dev.off()
###############

### e56 result plotting
###############
e56.umap <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e56-feat-umap.rds'); umap.pref <- 'fu'
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
exp5.samps <- colnames(e56.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK006')]
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p02-pwf.rds'); model.name <- 'k10-e56-dd-p02-pwf'
## classify
e56.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
e56.mvc <- e56.pred$aggregate[c(exp5.samps, exp6.samps)]
e56.percent.mat <- MakePercentMat(e56.pred)
class.colors <- ClusterColors(4); names(class.colors) <- colnames(e56.percent.mat)
## generate color vectors
class.thresh.col.vect <- ClassThreshColVect(e56.percent.mat, class.colors, class.thresh = 0.2)
class.entropy.vect <- ClassEntropyVect(e56.percent.mat)
class.alpha <- class.entropy.vect - min(class.entropy.vect); class.alpha <- class.alpha * (0.9 / max(class.alpha))
class.alpha <- 1 - class.alpha
## generate data frame
ak5.class <- as.vector(e56.mvc); names(ak5.class) <- names(e56.mvc); ak11.class[exp5.samps] <- 'AK005'
ak6.class <- as.vector(e56.mvc); names(ak6.class) <- names(e56.mvc); ak12.class[exp6.samps] <- 'AK006'
full.class.colors <- class.colors; full.class.colors[['AK005']] <- 'grey'; full.class.colors[['AK006']] <- 'grey'
plot.dat <- data.frame('UMAP1' = e56.umap$layout[,1], 'UMAP2' = e56.umap$layout[,2],
                       'class' = e56.mvc, 'colors' = class.thresh.col.vect, 'alpha' = I(class.alpha),
                       'entropy' = class.entropy.vect, 'sample' = substr(rownames(e56.umap$layout), 1, 5),
                       'ak5class' = ak5.class, 'ak6class' = ak6.class)
## sample specific data frames
ak5.alpha <- class.alpha[exp5.samps]; ak6.alpha <- class.alpha[exp6.samps]
ak5.dat <- data.frame('UMAP1' = e56.umap$layout[exp5.samps, 1], 'UMAP2' = e56.umap$layout[exp5.samps, 2],
                       'class' = e56.mvc[exp5.samps], 'alpha' = I(ak5.alpha),
                       'class.mix' = class.thresh.col.vect[exp5.samps])
ak6.dat <- data.frame('UMAP1' = e56.umap$layout[exp6.samps, 1], 'UMAP2' = e56.umap$layout[exp6.samps, 2],
                       'class' = e56.mvc[exp6.samps], 'alpha' = I(ak6.alpha),
                       'class.mix' = class.thresh.col.vect[exp6.samps])
ak5.dat <- cbind(ak5.dat, e56.percent.mat[exp5.samps,])
ak6.dat <- cbind(ak6.dat, e56.percent.mat[exp6.samps,])
x.min <- min(e56.umap$layout[,1]) - 0.5; x.max <- max(e56.umap$layout[,1]) + 0.5 
y.min <- min(e56.umap$layout[,2]) - 0.5; y.max <- max(e56.umap$layout[,2]) + 0.5
## separate mvc plots
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('Pooled') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
p2 <- ggplot(ak5.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('AK5') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
p3 <- ggplot(ak6.dat, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = class)) + 
  geom_density_2d(color = 'forestgreen') + ggtitle('AK6') +
  xlim(x.min, x.max) + ylim(y.min, y.max)
e56.mvc.sep <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e56-mvc-sep', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e56.mvc.sep, top = text_grob('AK005 / AK006 Classification', size = 16)))
dev.off()
## separate mvc plots w/ entropy alpha
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('Pooled') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen') 
p2 <- ggplot(ak5.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK5') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen') 
p3 <- ggplot(ak6.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK6') +
  geom_point(aes(color = class, alpha = alpha)) + geom_density_2d(color = 'forestgreen')
e56.mvc.sep.ent <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e56-mvc-sep-ent', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e56.mvc.sep.ent, top = text_grob('AK005 / AK006 Classification', size = 16)))
dev.off()
## separate mvc plots w/ color mixing
p1 <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('Pooled') +
  geom_point(aes(color = colors)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
p2 <- ggplot(ak5.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK5') +
  geom_point(aes(color = class.mix)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
p3 <- ggplot(ak6.dat, aes(x = UMAP1, y = UMAP2)) + xlim(x.min, x.max) + ylim(y.min, y.max) + ggtitle('AK6') +
  geom_point(aes(color = class.mix)) + geom_density_2d(color = 'forestgreen') +
  scale_color_identity()
e56.mvc.sep.col.mix <- ggarrange(plotlist = list('1' = p1, '2' = p2, '3' = p3), nrow = 1, ncol = 3)
jpeg(file = paste('plots/model-results/e56-mvc-sep-col-mix', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 1800, height = 700)
print(annotate_figure(e56.mvc.sep.col.mix, top = text_grob('AK005 / AK006 Classification', size = 16)))
dev.off()
## eight pannel plot
plot.list <- list()
for (i in names(class.colors)) {
  # AK11 plot
  plot.obj <- ggplot(ak5.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(ak5.dat[[i]])) + 
    ggtitle(paste('AK5: ', toupper(i), ' Classification', sep = '')) + xlim(x.min, x.max) + ylim(y.min, y.max) +
    geom_density_2d(color = 'forestgreen')
  plot.list[[paste('AK5', i, sep = '.')]] <- plot.obj
  # AK12 plot
  plot.obj <- ggplot(ak6.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(ak6.dat[[i]])) + 
    ggtitle(paste('AK6: ', toupper(i), ' Classification', sep = '')) + xlim(x.min, x.max) + ylim(y.min, y.max) +
    geom_density_2d(color = 'forestgreen')
  plot.list[[paste('AK6', i, sep = '.')]] <- plot.obj
}
plot.list <- plot.list[sort(names(plot.list))]
e56.eight.panel <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 4)
jpeg(file = paste('plots/model-results/e56-eight-panel', model.name, umap.pref, '.jpg', sep = '_'), 
     width = 2000, height = 1100)
print(annotate_figure(e56.eight.panel, top = text_grob('AK005 / AK006 Classification', family = 'Arial', size = 16)))
dev.off()
###############

### result plotting (original)
###############
e56.umap <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vumap.rds')
e1112.umap <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-1112-vumap.rds')
umap.pref <- 'u'
## load viper matrices
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
e56.lnorm <- e56.lnorm[, which(substr(colnames(e56.lnorm), 1, 5) %in% c('AK005', 'AK006'))]
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
## classify
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p02-pwf.rds')
model.name <- 'k10-dd-p02-pwf'
e56.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
k10rnv.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
## organize classifcation
e56nv.mvc <- e56.pred$aggregate
e56.percent.mat <- MakePercentMat(e56.pred)
e56nv.yi <- YIPredict(e56.percentMat, model.obj$model.yt)
e1112.mvc <- k10rnv.pred$aggregate; e1112.mvc <- e1112.mvc[which(substr(names(e1112.mvc), 1, 5) %in% c('AK011', 'AK012'))]
e1112.percent.mat <- MakePercentMat(k10rnv.pred)
e1112.yi <- YIPredict(e1112.percent.mat, model.obj$model.yt); e1112.yi <- e1112.yi[which(substr(names(e1112.yi), 1, 5) %in% c('AK011', 'AK012'))]
## get population names
exp5.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK006')]
exp11.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK012')]
## create plot matrices
e56.plot.dat <- data.frame('UMAP1' = e56.umap$layout[,1], 'UMAP2' = e56.umap$layout[,2],
                           'hsc' = e56.percent.mat[,1], 'mlp' = e56.percent.mat[,2], 'mpp' = e56.percent.mat[,3], 'prog' = e56.percent.mat[,4],
                           'class' = e56nv.mvc)
e1112.plot.dat <- data.frame('UMAP1' = e1112.umap$layout[,1], 'UMAP2' = e1112.umap$layout[,2],
                           'hsc' = e1112.percent.mat[c(exp11.samps, exp12.samps),1], 'mlp' = e1112.percent.mat[c(exp11.samps, exp12.samps),2], 
                           'mpp' = e1112.percent.mat[c(exp11.samps, exp12.samps),3], 'prog' = e1112.percent.mat[c(exp11.samps, exp12.samps),4],
                           'class' = e1112.mvc)
class.colors <- ClusterColors(4); names(class.colors) <- c('hsc', 'mlp', 'mpp', 'prog')
## vanilla umaps
jpeg(file = paste('plots/model-results/e56-nv_class', model.name, umap.pref, '.jpg', sep = '_'), width = 800, height = 900)
ggplot(e56.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK005 / AK006 Classification')
dev.off()
jpeg(file = paste('plots/model-results/e1112_class', model.name, umap.pref, '.jpg', sep = '_'), width = 800, height = 900)
ggplot(e1112.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK011 / AK012 Classification')
dev.off()
## e56 split MVC plot
e5.class <- as.character(e56nv.mvc); names(e5.class) <- names(e56nv.mvc); e5.class[exp6.samps] <- 'AK006'
e6.class <- as.character(e56nv.mvc); names(e6.class) <- names(e56nv.mvc); e6.class[exp5.samps] <- 'AK005'
e5.split.plot.dat <- data.frame('UMAP1' = e56.plot.dat$UMAP1, 'UMAP2' = e56.plot.dat$UMAP2,
                                 'class' = e5.class)
e6.split.plot.dat <- data.frame('UMAP1' = e56.plot.dat$UMAP1, 'UMAP2' = e56.plot.dat$UMAP2,
                                'class' = e6.class)
e5.cols <- class.colors; e5.cols[['AK006']] <- 'lightgrey'
e6.cols <- class.colors; e6.cols[['AK005']] <- 'lightgrey'
p1 <- ggplot(e5.split.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK005') +
  scale_color_manual(values = e5.cols)
p2 <- ggplot(e6.split.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK006 ') +
  scale_color_manual(values = e6.cols)
e56.split.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = paste('plots/model-results/e56-nv_split-class', model.name, umap.pref, '.jpg', sep = '_'), width = 1200, height = 700)
print(annotate_figure(e56.split.plot, top = text_grob('AK005 / AK006 Classification', family = 'Arial', size = 16)))
dev.off()
## e1112 split MVC plot
e11.class <- as.character(e1112.mvc); names(e11.class) <- names(e1112.mvc); e11.class[exp12.samps] <- 'AK012'
e12.class <- as.character(e1112.mvc); names(e12.class) <- names(e1112.mvc); e12.class[exp11.samps] <- 'AK011'
e11.split.plot.dat <- data.frame('UMAP1' = e1112.plot.dat$UMAP1, 'UMAP2' = e1112.plot.dat$UMAP2,
                                'class' = e11.class)
e12.split.plot.dat <- data.frame('UMAP1' = e1112.plot.dat$UMAP1, 'UMAP2' = e1112.plot.dat$UMAP2,
                                'class' = e12.class)
e11.cols <- class.colors; e11.cols[['AK012']] <- 'lightgrey'
e12.cols <- class.colors; e12.cols[['AK011']] <- 'lightgrey'
p1 <- ggplot(e11.split.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK011') +
  scale_color_manual(values = e11.cols)
p2 <- ggplot(e12.split.plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = class)) + ggtitle('AK012') +
  scale_color_manual(values = e12.cols)
e1112.split.plot <- ggarrange(plotlist = list('1' = p1, '2' = p2), nrow = 1, ncol = 2)
jpeg(file = paste('plots/model-results/e1112_split-class', model.name, umap.pref, '.jpg', sep = '_'), width = 1200, height = 700)
print(annotate_figure(e1112.split.plot, top = text_grob('AK011 / AK012 Classification', family = 'Arial', size = 16)))
dev.off()
## e56 eight pannel

## e1112 eight pannel
e1112.plot.title <- 'AK011 / AK012: Model Classification'
plot.list <- list()
for (i in names(class.colors)) { 
  ## AK011 plot
  plot.obj <- ggplot(e1112.plot.dat[exp11.samps,], aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(e1112.plot.dat[exp11.samps,][[i]])) + 
    ggtitle(paste('AK011: ', toupper(i), ' Classification', sep = ''))
  plot.list[[paste('AK011', i, sep = '.')]] <- plot.obj
  ## AK012 plot
  plot.obj <- ggplot(e1112.plot.dat[exp12.samps,], aes(UMAP1, UMAP2)) + geom_point(aes_string(color = i)) +
    scale_color_gradient2(low = 'darkgrey', high = class.colors[[i]], limits = range(e1112.plot.dat[exp12.samps,][[i]])) + 
    ggtitle(paste('AK012: ', toupper(i), ' Classification', sep = ''))
  plot.list[[paste('AK012', i, sep = '.')]] <- plot.obj
}
plot.list <- plot.list[sort(names(plot.list))]
e1112.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 4)
jpeg(file = paste('plots/model-results/e1112_fuzzy-class', model.name, umap.pref, '.jpg', sep = '_'), width = 2000, height = 1100)
print(annotate_figure(e1112.plot, top = text_grob(e1112.plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### maker VIPER stream objects
###############
exp.vip <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
e56.vip <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p02.rds')
## get population names
exp5.samps <- colnames(e56.vip)[which(substr(colnames(e56.vip), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.vip)[which(substr(colnames(e56.vip), 1, 5) == 'AK006')]
exp11.samps <- colnames(exp.vip)[which(substr(colnames(exp.vip), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.vip)[which(substr(colnames(exp.vip), 1, 5) == 'AK012')]
## save stream data matrices
write.table(e56.vip[, exp5.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e5-vip_stream.tsv.gz'))
write.table(e56.vip[, exp6.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e6-vip_stream.tsv.gz'))
write.table(exp.vip[, exp11.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e11-vip_stream.tsv.gz'))
write.table(exp.vip[, exp12.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e12-vip_stream.tsv.gz'))
###############

### reference stream objects
###############
ref.filt <- readRDS('data/ref-samps/ref-pops_filt.rds')
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p02.rds')
## write matrices
write.table(ref.vip[, names(train.vect)], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/ref-vip_stream.tsv.gz'))
write.table(ref.filt[, names(train.vect)], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/ref_stream.tsv.gz'))
## write metadata
class.colors <- ClusterColors(4)
train.vect <- as.factor(train.vect)
meta.df <- data.frame('label' = train.vect, 'label_color' = mapvalues(train.vect, levels(train.vect), class.colors))
write.table(meta.df, sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/ref_dd-p02-samps.tsv.gz'))
###############

### make stream objects
###############
exp.filt <- readRDS('data/exp-samps/exp-pops_filt.rds')
ref.filt <- readRDS('data/ref-samps/ref-pops_filt.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p02.rds')
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p02-pwf.rds')
## load viper and classify
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
e56.lnorm <- e56.lnorm[, which(substr(colnames(e56.lnorm), 1, 5) %in% c('AK005', 'AK006'))]
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
e56.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
k10rnv.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
## get population names
exp5.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK005')]
exp6.samps <- colnames(e56.lnorm)[which(substr(colnames(e56.lnorm), 1, 5) == 'AK006')]
exp11.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK011')]
exp12.samps <- colnames(exp.lnorm)[which(substr(colnames(exp.lnorm), 1, 5) == 'AK012')]
## save stream data matrices
write.table(exp.filt[, exp5.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e5_stream.tsv.gz'))
write.table(exp.filt[, exp6.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e6_stream.tsv.gz'))
write.table(exp.filt[, exp11.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e11_stream.tsv.gz'))
write.table(exp.filt[, exp12.samps], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
            file = gzfile('data/stream-data/e12_stream.tsv.gz'))
## save stream metadata objects
class.vect <- c(as.character(e56.pred$aggregate[exp5.samps]), as.character(e56.pred$aggregate[exp6.samps]),
                as.character(k10rnv.pred$aggregate[exp11.samps]), as.character(k10rnv.pred$aggregate[exp12.samps]))
names(class.vect) <- c(exp5.samps, exp6.samps, exp11.samps, exp12.samps); class.vect <- as.factor(class.vect)
class.colors <- ClusterColors(4)
meta.df <- data.frame('label' = class.vect, 'label_color' = mapvalues(class.vect, levels(class.vect), class.colors))
rownames(meta.df) <- names(class.vect)
write.table(meta.df[exp5.samps,], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
          file = gzfile('data/stream-data/e5_e56nv-k10rnv-dd-p02-pwf-mvc.tsv.gz'))
write.table(meta.df[exp6.samps,], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
          file = gzfile('data/stream-data/e6_e56nv-k10rnv-dd-p02-pwf-mvc.tsv.gz'))
write.table(meta.df[exp11.samps,], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
          file = gzfile('data/stream-data/e11_k10rnv-dd-p02-pwf-mvc.tsv.gz'))
write.table(meta.df[exp12.samps,], sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE,
          file = gzfile('data/stream-data/e12_k10rnv-dd-p02-pwf-mvc.tsv.gz'))
###############

### progenitor population mrs + markers
###############
vip.mat <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
vip.umap <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vumap.rds')
samp.labels <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p02.rds')
## convert to gene name
dim(vip.mat)
vip.mat <- Ensemble2GeneName(vip.mat)
dim(vip.mat)
## identify MRS for progenitor population
prog.samps <- names(samp.labels)[which(samp.labels == 'prog')]; n1 <- length(prog.samps)
ref.samps <- names(samp.labels)[which(samp.labels != 'prog')]; n2 <- length(ref.samps)
w.stats <- apply(vip.mat, 1, function(x) {
  return(wilcox.test(x[prog.samps], x[ref.samps], alternative = 'greater')$statistic)
})
w.mean <- (n1 * n2) / 2
w.sd <- sqrt(((n1 * n2) * (n1 + n2 + 1)) / 12)
p.vals <- unlist(lapply(w.stats, function(x) { pnorm(((x - w.mean) / w.sd), lower.tail = FALSE, log.p = TRUE) }))
p.vals <- sort(p.vals)
write.table(as.data.frame(p.vals[1:100]), file = 'model-results/k10rnv_dd-p02_ref-prog-mrs.tsv',
            sep = '\t', col.names = FALSE, row.names = TRUE, quote = FALSE)
## set up poplation labels
plot.labels <- rep('NA', nrow(vip.umap$layout)); names(plot.labels) <- rownames(vip.umap$layout)
for (u in unique(samp.labels)) {
  plot.labels[names(samp.labels)[which(samp.labels == u)]] <- u
}
## marker plot w/ data frame
markers <- c('GATA2', 'FOXM1', 'MYBL2', 'GATA1', 'KLF7')
markers <- intersect(rownames(vip.mat), markers)
plot.dat <- data.frame('UMAP1' = vip.umap$layout[,1], 'UMAP2' = vip.umap$layout[,2],
                       'label' = plot.labels)
plot.dat <- cbind(plot.dat, as.data.frame(t(vip.mat[markers,])))
class.cols <- c(ClusterColors(4), 'grey'); names(class.cols) <- c('hsc', 'mpp', 'mlp', 'prog', 'NA')
## generate plots
plot.list <- list()
plot.list[['pops']] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = label)) + ggtitle('Training Populations') +
  scale_color_manual(values = class.cols) + geom_density_2d(color = 'forestgreen')
for (m in markers) {
  plot.obj <- ggplot(plot.dat, aes_string(x = 'UMAP1', y = 'UMAP2', color = m)) + geom_point() + ggtitle(m) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', na.value = 'lightgray', midpoint = mean(plot.dat[,m])) 
  plot.list[[m]] <- plot.obj
}
plot.title <- 'Training Populations: HSC Markers'
m.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = 3)
jpeg(file = 'plots/model-results/k10rnv-p02_ref-train-samps_markers.jpg', width = 1500, height = 1100)
print(annotate_figure(m.plot, top = text_grob(plot.title, family = 'Arial', size = 16)))
dev.off()
###############

### AK11 v AK12 gene expression
###############
exp.cpm <- readRDS('data/exp-samps/exp-pops_cpm.rds')
exp.cpm <- Ensemble2GeneName(exp.cpm)
exp.filt <- readRDS('data/exp-samps/exp-pops_filt.rds')
exp.filt <- Ensemble2GeneName(exp.filt)
samp.names <- substr(colnames(exp.cpm), 1, 5)
ak11.samps <- colnames(exp.cpm)[which(samp.names == 'AK011')]; n1 <- length(ak11.samps)
ak12.samps <- colnames(exp.cpm)[which(samp.names == 'AK012')]; n2 <- length(ak12.samps)
## get mean and variance
w.mean <- (n1 * n2) / 2
w.sd <- sqrt(((n1 * n2) * (n1 + n2 + 1)) / 12)
## get ak11 p-values
w.stats <- apply(exp.cpm, 1, function(x) {
  return(wilcox.test(x[ak11.samps], x[ak12.samps], alternative = 'greater')$statistic)
})
p.vals <- sapply(w.stats, function(x) { pnorm(((x - w.mean) / w.sd), lower.tail = FALSE, log.p = TRUE) })
p.vals <- sort(p.vals)
write.table(as.data.frame(p.vals[1:100]), file = 'data/exp-samps/ak11_dif-gexp.tsv', 
            col.names = FALSE, row.names = TRUE, quote = FALSE, sep = '\t')
## get ak12 p-values
w.stats <- apply(exp.cpm, 1, function(x) {
  return(wilcox.test(x[ak12.samps], x[ak11.samps], alternative = 'greater')$statistic)
})
p.vals <- sapply(w.stats, function(x) { pnorm(((x - w.mean) / w.sd), lower.tail = FALSE, log.p = TRUE) })
p.vals <- sort(p.vals)
write.table(as.data.frame(p.vals[1:100]), file = 'data/exp-samps/ak12_dif-gexp.tsv', 
            col.names = FALSE, row.names = TRUE, quote = FALSE, sep = '\t')
## remove genes with zero counts
ak11.genes <- rownames(exp.filt)[which(rowSums(exp.filt[,ak11.samps]) > 0)]
ak12.genes <- rownames(exp.filt)[which(rowSums(exp.filt[,ak12.samps]) > 0)]
exp.filt <- exp.filt[intersect(ak11.genes, ak12.genes),]
## bootstrap log fold changes
boot.num <- 100
lfc.mat <- matrix(0L, nrow = nrow(exp.filt), ncol = boot.num)
rownames(lfc.mat) <- rownames(exp.filt); colnames(lfc.mat) <- 1:boot.num
for (b in 1:boot.num) {
  print(b)
  # generate bootstrap
  ak11.boot <- rowSums(exp.filt[, sample(ak11.samps, length(ak11.samps), replace = TRUE)]) + 1
  ak11.boot <- (ak11.boot / sum(ak11.boot)) * 1e6
  ak12.boot <- rowSums(exp.filt[, sample(ak12.samps, length(ak12.samps), replace = TRUE)]) + 1
  ak12.boot <- (ak12.boot / sum(ak12.boot)) * 1e6
  # calculate and store
  lfc.mat[,b] <- log(ak12.boot / ak11.boot, base = 2)
}
## bootstrap means / sd
lfc.means <- rowMeans(lfc.mat)
lfc.sd <- apply(lfc.mat, 1, sd)
lfc.df <- data.frame('Mean' = lfc.means, 'SD' = lfc.sd)
rownames(lfc.df) <- rownames(lfc.mat)
lfc.df <- lfc.df[order(lfc.df$Mean, decreasing = TRUE),]
write.table(lfc.df, file = 'data/exp-samps/ak12v11_gexp-lfc.tsv', 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
saveRDS(lfc.df, file = 'data/exp-samps/ak12v11_gexp-lfc.rds')
###############

### LFC + p-val analysis w/ gExp
###############
ref.cpm <- readRDS('data/ref-samps/ref-pops_cpm.rds')
exp.cpm <- readRDS('data/exp-samps/exp-pops_cpm.rds')
ak1112.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p01-pwf_e1112-mvc.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p01.rds')
## get exp names
hsc.names <- names(ak1112.mvc)[which(ak1112.mvc == 'hsc')]
prog.names <- names(ak1112.mvc)[which(ak1112.mvc == 'prog')]
ak11.names <- colnames(exp.cpm)[which(substr(colnames(exp.cpm), 1, 5) == 'AK011')]
ak12.names <- colnames(exp.cpm)[which(substr(colnames(exp.cpm), 1, 5) == 'AK012')]
## training hsc v prog
a.mat <- ref.cpm[,names(train.vect)[which(train.vect == 'hsc')]]
b.mat <- ref.cpm[,names(train.vect)[which(train.vect == 'prog')]]
s.genes <- intersect(rownames(a.mat)[which(rowSums(a.mat) > 0)], rownames(b.mat)[which(rowSums(b.mat) > 0)])
gene.lfc <- ISBLogFoldChange(a.mat[s.genes,], b.mat[s.genes,], boot.num = 100)
gene.pVal <- c(); shared.genes <- intersect(rownames(a.mat), rownames(b.mat))
for (g in s.genes) {
  gene.pVal <- c(gene.pVal, wilcox.test(a.mat[g,], b.mat[g,], alternative = 'two.sided', exact = TRUE)$p.val)
}
names(gene.pVal) <- s.genes; gene.pVal <- gene.pVal[rownames(gene.lfc)]
gene.lfc[['pVal']] <- gene.pVal
saveRDS(gene.lfc, file = 'data/gExp-analysis/ref_dd-p01_prog-v-hsc_lfc.rds')
## AK12 v. AK11
a.mat <- exp.cpm[,ak11.names]
b.mat <- exp.cpm[,ak12.names]
s.genes <- intersect(rownames(a.mat)[which(rowSums(a.mat) > 0)], rownames(b.mat)[which(rowSums(b.mat) > 0)])
gene.lfc <- ISBLogFoldChange(a.mat[s.genes,], b.mat[s.genes,], boot.num = 100)
gene.pVal <- c(); shared.genes <- intersect(rownames(a.mat), rownames(b.mat))
for (g in s.genes) {
  gene.pVal <- c(gene.pVal, wilcox.test(a.mat[g,], b.mat[g,], alternative = 'two.sided', exact = TRUE)$p.val)
}
names(gene.pVal) <- s.genes; gene.pVal <- gene.pVal[rownames(gene.lfc)]
gene.lfc[['pVal']] <- gene.pVal
saveRDS(gene.lfc, file = 'data/gExp-analysis/exp_ak12-v-ak11_lfc.rds')
## AK1112 hsc v. AK1112 prog
a.mat <- exp.cpm[,hsc.names]
b.mat <- exp.cpm[,prog.names]
s.genes <- intersect(rownames(a.mat)[which(rowSums(a.mat) > 0)], rownames(b.mat)[which(rowSums(b.mat) > 0)])
gene.lfc <- ISBLogFoldChange(a.mat[s.genes,], b.mat[s.genes,], boot.num = 100)
gene.pVal <- c(); shared.genes <- intersect(rownames(a.mat), rownames(b.mat))
for (g in s.genes) {
  gene.pVal <- c(gene.pVal, wilcox.test(a.mat[g,], b.mat[g,], alternative = 'two.sided', exact = TRUE)$p.val)
}
names(gene.pVal) <- s.genes; gene.pVal <- gene.pVal[rownames(gene.lfc)]
gene.lfc[['pVal']] <- gene.pVal
saveRDS(gene.lfc, file = 'data/gExp-analysis/exp_dd-p01-pwf-mvc_prog-v-hsc_lfc.rds')
## AK12 hsc v. AK11 prog
a.mat <- exp.cpm[,intersect(hsc.names, ak12.names)]
b.mat <- exp.cpm[,intersect(prog.names, ak11.names)]
s.genes <- intersect(rownames(a.mat)[which(rowSums(a.mat) > 0)], rownames(b.mat)[which(rowSums(b.mat) > 0)])
gene.lfc <- ISBLogFoldChange(a.mat[s.genes,], b.mat[s.genes,], boot.num = 100)
gene.pVal <- c(); shared.genes <- intersect(rownames(a.mat), rownames(b.mat))
for (g in s.genes) {
  gene.pVal <- c(gene.pVal, wilcox.test(a.mat[g,], b.mat[g,], alternative = 'two.sided', exact = TRUE)$p.val)
}
names(gene.pVal) <- s.genes; gene.pVal <- gene.pVal[rownames(gene.lfc)]
gene.lfc[['pVal']] <- gene.pVal
saveRDS(gene.lfc, file = 'data/gExp-analysis/exp_dd-p01-pwf-mvc_ak11prog-v-ak12hsc_lfc.rds')
###############

### Reference population LFC
###############
ref.filt <- readRDS('data/ref-samps/ref-pops_filt.rds')
ref.filt <- Ensemble2GeneName(ref.filt)
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p02.rds')
## subset hsc and prog
hsc.mat <- ref.filt[, which(train.vect == 'hsc')]
prog.mat <- ref.filt[, which(train.vect == 'prog')]
lfc.df <- ISBLogFoldChange(hsc.mat, prog.mat)
write.table(lfc.df, file = 'data/ref-samps/k10rnv-dd-p02_hsc-v-prog_lfc.tsv', 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
saveRDS(lfc.df, file = 'data/ref-samps/k10rnv-dd-p02_hsc-v-prog_lfc.rds')
###############

## class unit circle plot
###############
model.name <- 'k10-e56-dd-p02-pwf'
e56.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e56-mvc.rds')
e56.percent.mat <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e56-percent-mat.rds')
e1112.mvc <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e1112-mvc.rds')
e1112.percent.mat <- readRDS('model-results/classification-objects/k10-e56-dd-p02-pwf_e1112-percent-mat.rds')
## calculate entropy
e56.entropy <- ClassEntropyVect(e56.percent.mat)
e1112.entropy <- ClassEntropyVect(e1112.percent.mat)

## set params
ent.vect <- e56.entropy
mvc.vect <- e56.mvc
per.mat <- e56.percent.mat

# plot parameters
plot.rad <- 10; min.rad <- 1; max.ent <- log(4, base = 2)
class.angles <- c(pi, 3*pi, 5*pi, 7*pi) / 4; names(class.angles) <- c('hsc', 'mlp', 'mpp', 'prog')
# adjust entropy and calculate angles
ent.vect <- (-1) * (ent.vect - max(ent.vect)) # invert so maximal entropy is at the origin
ent.vect <- (plot.rad / max.ent) * ent.vect # rescale so that the maximum possible entropy would lie on the circle
#samp.angles <- rowSums(t(t(per.mat) * class.angles))
squared.per.mat <- per.mat**2; squared.per.mat <- squared.per.mat / rowSums(squared.per.mat); samp.angles <- rowSums(t(t(squared.per.mat) * class.angles))
# carteisan transformation
x.coords <- ent.vect * cos(samp.angles)
y.coords <- ent.vect * sin(samp.angles)
plot.dat <- data.frame('x.coord' = x.coords, 'y.coord' = y.coords, 'class' = mvc.vect)

## plotting
gg_circle <- function(r, xc, yc, color="black", fill='NA', ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}
circle.int <- plot.rad / sqrt(2)
class.cols <- ClusterColors(4); names(class.cols) <- names(class.angles)
text.buffer <- 1; anno.size <- 8

plot.title <- 'AK005'
plot.samps <- which(substr(rownames(plot.dat), 1, 5) == plot.title)
jpeg(file = paste('plots/model-results/', tolower(plot.title), '-circle-plot-2_', model.name, '.jpg', sep = ''),  width = 1100, height = 1050)
ggplot(plot.dat[plot.samps,], aes(x.coord, y.coord)) + gg_circle(r = 10, xc = 0, yc = 0, fill = 'white') + geom_point(aes(color = class)) + 
  geom_segment(x = -circle.int, y = -circle.int, xend = circle.int, yend = circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
  geom_segment(x = -circle.int, y = circle.int, xend = circle.int, yend = -circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
  xlim(c(-10,10)) + ylim(c(-10,10)) + ggtitle(plot.title) + geom_density_2d(color = 'forestgreen') +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), text = element_text(size = 16), ) +
  annotate(geom="text", x = circle.int + text.buffer, y = circle.int + text.buffer, label="HSC", color = class.cols[['hsc']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = -circle.int - text.buffer, y = circle.int + text.buffer, label="MLP", color = class.cols[['mlp']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = -circle.int - text.buffer, y =  -circle.int - text.buffer, label="MPP", color = class.cols[['mpp']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = circle.int + text.buffer, y = -circle.int - text.buffer, label="PROG", color = class.cols[['prog']], fontface = 'bold', size = anno.size)
dev.off()

plot.title <- 'AK006'
plot.samps <- which(substr(rownames(plot.dat), 1, 5) == plot.title)
jpeg(file = paste('plots/model-results/', tolower(plot.title), '-circle-plot-2_', model.name, '.jpg', sep = ''),  width = 1100, height = 1050)
ggplot(plot.dat[plot.samps,], aes(x.coord, y.coord)) + gg_circle(r = 10, xc = 0, yc = 0, fill = 'white') + geom_point(aes(color = class)) + 
  geom_segment(x = -circle.int, y = -circle.int, xend = circle.int, yend = circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
  geom_segment(x = -circle.int, y = circle.int, xend = circle.int, yend = -circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
  xlim(c(-10,10)) + ylim(c(-10,10)) + ggtitle(plot.title) + geom_density_2d(color = 'forestgreen') +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), text = element_text(size = 16)) +
  annotate(geom="text", x = circle.int + text.buffer, y = circle.int + text.buffer, label="HSC", color = class.cols[['hsc']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = -circle.int - text.buffer, y = circle.int + text.buffer, label="MLP", color = class.cols[['mlp']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = -circle.int - text.buffer, y =  -circle.int - text.buffer, label="MPP", color = class.cols[['mpp']], fontface = 'bold', size = anno.size) +
  annotate(geom="text", x = circle.int + text.buffer, y = -circle.int - text.buffer, label="PROG", color = class.cols[['prog']], fontface = 'bold', size = anno.size)
dev.off()
###############

## identify prog specific features
###############
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p01.rds')
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p01-pwf.rds')
model.feats <- model.obj$model.feats
## get class specific p-values for each feature
dat.labels <- as.factor(train.vect)
dat.classes <- levels(dat.labels)
pVal.mat <- apply(ref.vip, 1, function(x) {
  p.mat <- pairwise.wilcox.test(x, dat.labels)$p.value
  p.vals <- list()
  for (cl in dat.classes) {
    class.p <- c(p.mat[, intersect(colnames(p.mat), cl)], p.mat[ intersect(rownames(p.mat), cl) ,])
    class.p <- class.p[ !is.na(class.p) ]
    p.vals[[cl]] <- min(class.p)
  }
  return(unlist(p.vals))
})
## get genes for which prog had the lowest p.value
sorted.mat <- apply(pVal.mat, 1, function(x) {
  return(names(sort(x)))
})
## find and save prog genes
prog.genes <- intersect(sorted.mat[1:20,4], model.feats)
write.table(prog.genes, file = 'model-results/dd-p01-pwf_prog-genes.tsv', sep = '\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE)
###############

## reference stouffer MRs
###############
ref.vip <- readRDS('data/viper-mats/k10rnv/ref-pops_k10rn-vip.rds')
train.vect <- readRDS('model-results/refined-labels/k10rnv_dif-density-labels_p01.rds')
## get MRS for each group
mr.list <- list()
thresh.sig <- 10
for (p in unique(train.vect)) {
  p.samps <- names(train.vect)[which(train.vect == p)]
  p.mat <- ref.vip[, p.samps]
  p.mrs <- StoufferMRs(p.mat)
  mr.list[[p]] <- sort(p.mrs, decreasing = TRUE)[1:300]
}
saveRDS(mr.list, file = 'model-results/k10rnv_dd-p01_mrs.rds')
###############

## model feature gene names
###############
model.obj <- readRDS('model-results/model-objs/k10rnv_dd-p01-pwf.rds')
model.feats <- model.obj$model.feats
## convert to gene name
num.feats <- length(model.feats)
dummy.mat <- matrix(1:(2*num.feats), ncol = 2); rownames(dummy.mat) <- model.feats
dummy.mat <- GeneNameConvert(dummy.mat, 'hum', 'ensg', 'gn')
model.feats.gn <- data.frame('ENSG' = model.feats, 'GN' = rownames(dummy.mat))
## write table
write.table(model.feats.gn, file = 'model-results/k10rnv_dd-p01-pwf_feats.tsv', sep = '\t',
            quote = FALSE, col.names = TRUE, row.names = FALSE)
###############

## atac TF activity
###############
vip.mat <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
tf.mat <- read.csv('atac/atac-tf-short-list.csv', stringsAsFactors = FALSE)
tf.list <- unlist(tf.mat, use.names = FALSE)
tf.list <- tf.list[which(tf.list != '')]
## stouffer integrate
vip.mat.gn <- GeneNameConvert(vip.mat, species = 'hum', 'ensg', 'gn')
ak12.mat <- vip.mat.gn[, which(substr(colnames(vip.mat.gn), 1, 5) == 'AK012')]
## stouffer integrate 
sint.vect <- rowSums(vip.mat.gn)
sint.vect <- sint.vect / sqrt(ncol(vip.mat.gn))
sint.vect <- sort(sint.vect, decreasing = TRUE)
## get indices of tf.list
tf.inds <- match(tf.list, names(sint.vect))
names(tf.inds) <- tf.list
## save objects
saveRDS(tf.inds, file = 'atac/atac_tf-stouffer-inds.rds')
saveRDS(sint.vect, file = 'atac/k10rnv_ak12-stouffer-vect.rds')
###############



