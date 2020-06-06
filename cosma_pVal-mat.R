setwd('C://Users/lvlah/linux/ac_lab/cosma/')
library(randomForest)
library(plyr)
library(pROC)
library(MASS)
source('cosma_functions.R')


## load e56 matrices and subset for memory space
e56.rnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-rnorm-vip.rds')
e56.rnorm <- e56.rnorm[, which(substr(colnames(e56.rnorm), 1, 5) %in% c('AK005', 'AK006'))]
e56.lnorm <- readRDS('data/viper-mats/e56nv/e56_ep56n-lnorm-vip.rds')
e56.lnorm <- e56.lnorsavem[, which(substr(colnames(e56.lnorm), 1, 5) %in% c('AK005', 'AK006'))]
e56.rnorm.labels <- substr(colnames(e56.rnorm), 1, 5)
e56.lnorm.labels <- substr(colnames(e56.lnorm), 1, 5)
## load network-specific rnorm / lnorm
exp.rnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-rnorm-vip.rds')
exp.lnorm <- readRDS('data/viper-mats/k10rnv/exp-pops_k10rn-lnorm-vip.rds')
rnorm.labels <- substr(colnames(exp.rnorm), 1, 5)
lnorm.labels <- substr(colnames(exp.lnorm), 1, 5)

## get list of model objects for k10rnv
model.dir <- 'model-results/model-objs'
model.objs <- list.files(model.dir)
model.objs <- model.objs[!grepl('bcnv', model.objs)]
## create matrix of results
pVal.mat <- matrix(0L, nrow = length(model.objs), ncol = 12)
rownames(pVal.mat) <- unlist(lapply(model.objs, function(x) {strsplit(x, '\\.')[[1]][1]}))
colnames(pVal.mat) <- c('e56_lnorm_mvc', 'e56_lnorm_yi', 'e56_rnowm_mvc', 'e56_rnowm_yi',
                        'e56nv_lnorm_mvc', 'e56nv_lnorm_yi', 'e56nv_rnorm_mvc', 'e56nv_rnorm_yi',
                        'e1112_lnorm_mvc', 'e1112_lnmorm_yi', 'e1112_rnorm_mvc', 'e1112_rnorm_yi')

## loop through model list
for (m in model.objs) {
  # get name of model + full path
  print(m)
  model.name <- strsplit(m, '\\.')[[1]][1]
  model.path <- paste(model.dir, m, sep = '/')
  model.obj <- readRDS(model.path)
  # e56 lnorm
  exp.pred <- predict(model.obj$model, t(e56.lnorm[model.obj$model.feats,]), predict.all = TRUE)
  percent.mat <- MakePercentMat(exp.pred)
  yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
  if ('hsc' %in% exp.pred$aggregate) {
    mvc.hsc <- exp.pred$aggregate; mvc.hsc <- mapvalues(mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e56nv.lnorm.mvc.pval <- fisher.test(mvc.hsc, e56.lnorm.labels, alternative = 'less')$p.value
  } else {
    e56nv.lnorm.mvc.pval <- 1
  }
  if ('hsc' %in% yi.pred) {
    yi.hsc <- yi.pred; yi.hsc <- mapvalues(yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e56nv.lnorm.yi.pval <- fisher.test(yi.hsc, e56.lnorm.labels, alternative = 'less')$p.value
  } else {
    e56nv.lnorm.yi.pval <- 1
  }
  # e56 rnorm
  exp.pred <- predict(model.obj$model, t(e56.rnorm[model.obj$model.feats,]), predict.all = TRUE)
  percent.mat <- MakePercentMat(exp.pred)
  yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
  if ('hsc' %in% exp.pred$aggregate) {
    mvc.hsc <- exp.pred$aggregate; mvc.hsc <- mapvalues(mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e56nv.rnorm.mvc.pval <- fisher.test(mvc.hsc, e56.rnorm.labels, alternative = 'less')$p.value
  } else {
    e56nv.rnorm.mvc.pval <- 1
  }
  if ('hsc' %in% yi.pred) {
    yi.hsc <- yi.pred; yi.hsc <- mapvalues(yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e56nv.rnorm.yi.pval <- fisher.test(yi.hsc, e56.rnorm.labels, alternative = 'less')$p.value
  } else {
    e56nv.rnorm.yi.pval <- 1
  }
  ## lnorm test
  exp.pred <- predict(model.obj$model, t(exp.lnorm[model.obj$model.feats,]), predict.all = TRUE)
  percent.mat <- MakePercentMat(exp.pred)
  yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
  # e56
  e56.mvc.hsc <- exp.pred$aggregate; e56.mvc.hsc <- e56.mvc.hsc[which(substr(names(e56.mvc.hsc), 1, 5) %in% c('AK005', 'AK006'))]
  if ('hsc' %in% e56.mvc.hsc) {
    e56.mvc.hsc <- mapvalues(e56.mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e56.lnorm.mvc.pval <- fisher.test(e56.mvc.hsc, as.factor(substr(names(e56.mvc.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e56.lnorm.mvc.pval <- 1
  }
  e56.yi.hsc <- yi.pred; e56.yi.hsc <- e56.yi.hsc[which(substr(names(e56.yi.hsc), 1, 5) %in% c('AK005', 'AK006'))]
  if ('hsc' %in% e56.yi.hsc) {
    e56.yi.hsc <- mapvalues(e56.yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e56.lnorm.yi.pval <- fisher.test(e56.yi.hsc, as.factor(substr(names(e56.yi.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e56.lnorm.yi.pval <- 1
  }
  # e1112
  e1112.mvc.hsc <- exp.pred$aggregate; e1112.mvc.hsc <- e1112.mvc.hsc[which(substr(names(e1112.mvc.hsc), 1, 5) %in% c('AK011', 'AK012'))]
  if ('hsc' %in% e1112.mvc.hsc) {
    e1112.mvc.hsc <- mapvalues(e1112.mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e1112.lnorm.mvc.pval <- fisher.test(e1112.mvc.hsc, as.factor(substr(names(e1112.mvc.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e1112.lnorm.mvc.pval <- 1
  }
  e1112.yi.hsc <- yi.pred; e1112.yi.hsc <- e1112.yi.hsc[which(substr(names(e1112.yi.hsc), 1, 5) %in% c('AK011', 'AK012'))]
  if ('hsc' %in% e1112.yi.hsc) {
    e1112.yi.hsc <- mapvalues(e1112.yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e1112.lnorm.yi.pval <- fisher.test(e1112.yi.hsc, as.factor(substr(names(e1112.yi.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e1112.lnorm.yi.pval <- 1
  }
  ## rnorm test
  exp.pred <- predict(model.obj$model, t(exp.rnorm[model.obj$model.feats,]), predict.all = TRUE)
  percent.mat <- MakePercentMat(exp.pred)
  yi.pred <- YIPredict(percent.mat, model.obj$model.yt)
  # e56
  e56.mvc.hsc <- exp.pred$aggregate; e56.mvc.hsc <- e56.mvc.hsc[which(substr(names(e56.mvc.hsc), 1, 5) %in% c('AK005', 'AK006'))]
  if ('hsc' %in% e56.mvc.hsc) {
    e56.mvc.hsc <- mapvalues(e56.mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e56.rnorm.mvc.pval <- fisher.test(e56.mvc.hsc, as.factor(substr(names(e56.mvc.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e56.rnorm.mvc.pval <- 1
  }
  e56.yi.hsc <- yi.pred; e56.yi.hsc <- e56.yi.hsc[which(substr(names(e56.yi.hsc), 1, 5) %in% c('AK005', 'AK006'))]
  if ('hsc' %in% e56.yi.hsc) {
    e56.yi.hsc <- mapvalues(e56.yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e56.rnorm.yi.pval <- fisher.test(e56.yi.hsc, as.factor(substr(names(e56.yi.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e56.rnorm.yi.pval <- 1
  }
  # e1112
  e1112.mvc.hsc <- exp.pred$aggregate; e1112.mvc.hsc <- e1112.mvc.hsc[which(substr(names(e1112.mvc.hsc), 1, 5) %in% c('AK011', 'AK012'))]
  if ('hsc' %in% e1112.mvc.hsc) {
    e1112.mvc.hsc <- mapvalues(e1112.mvc.hsc, c('mpp', 'mlp', 'prog'), c('other', 'other', 'other'))
    e1112.rnorm.mvc.pval <- fisher.test(e1112.mvc.hsc, as.factor(substr(names(e1112.mvc.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e1112.rnorm.mvc.pval <- 1
  }
  e1112.yi.hsc <- yi.pred; e1112.yi.hsc <- e1112.yi.hsc[which(substr(names(e1112.yi.hsc), 1, 5) %in% c('AK011', 'AK012'))]
  if ('hsc' %in% e1112.yi.hsc) {
    e1112.yi.hsc <- mapvalues(e1112.yi.hsc, c('mpp', 'mlp', 'prog', 'unknown'), c('other', 'other', 'other', 'other'))
    e1112.rnorm.yi.pval <- fisher.test(e1112.yi.hsc, as.factor(substr(names(e1112.yi.hsc), 1, 5)), alternative = 'less')$p.value
  } else {
    e1112.rnorm.yi.pval <- 1
  }
  ## assemble p-val vector
  pVal.vect <- c(e56.lnorm.mvc.pval, e56.lnorm.yi.pval, e56.rnorm.mvc.pval, e56.rnorm.yi.pval,
                 e56nv.lnorm.mvc.pval, e56nv.lnorm.yi.pval, e56nv.rnorm.mvc.pval, e56nv.rnorm.yi.pval,
                 e1112.lnorm.mvc.pval, e1112.lnorm.yi.pval, e1112.rnorm.mvc.pval, e1112.rnorm.yi.pval)
  pVal.mat[model.name,] <- pVal.vect
}

write.table(pVal.mat, file = 'model-results/pVal-mat_2.csv', sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(pVal.mat, file = 'model-results/pVal-mat_2.rds')




