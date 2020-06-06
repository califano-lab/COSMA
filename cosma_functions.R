#' Performs feature selection on the given matrix and labels using Kruskal Wallis
#' 
#' @param dat.mat Matrix of data (features X samples)
#' @param dat.labels Factor of class labels
#' @param min.feats Minimum number of features. Default of 2.
#' @param num.iter Number of feature sets to return. Default of 49
#' @return Lists of lists; each sub-list is a candidate feature set.
kwFeats <- function(dat.mat, dat.labels, min.feats = 2, num.iter = 50) {
  # perform KW test
  dat.labels <- as.factor(dat.labels)
  kw.feats <- apply(dat.mat, 1, function(x) {
    return(kruskal.test(x, dat.labels)$p.value)
  })
  kw.feats <- sort(kw.feats)
  # generate list of feature sets
  feat.sets <- list()
  for (i in 1:num.iter) {
    feat.sets[[i]] <- head(names(kw.feats), (min.feats - 1 + i))
  }
  return(feat.sets)
}

#' Performs feature selection using a pairwise wilcox test.
#' For each gene, a per-class p-value is selected by taking the minimum p-value from all pairwise tests.
#' 
#' @param dat.mat Matrix of data (features X samples)
#' @param dat.labels Vector or factor of class labels
#' @param num.iter Number of candiate feature sets to return. Default of 25.
#' @return Lists of lists; each sub-list is a candidate feature set.
pwFeats <- function(dat.mat, dat.labels, num.iter = 25) {
  dat.labels <- as.factor(dat.labels)
  dat.classes <- levels(dat.labels)
  # calculate p.values for each gene
  pVal.mat <- apply(dat.mat, 1, function(x) {
    p.mat <- pairwise.wilcox.test(x, dat.labels)$p.value
    p.vals <- list()
    for (cl in dat.classes) {
      class.p <- c(p.mat[, intersect(colnames(p.mat), cl)], p.mat[ intersect(rownames(p.mat), cl) ,])
      class.p <- class.p[ !is.na(class.p) ]
      p.vals[[cl]] <- min(class.p)
    }
    return(unlist(p.vals))
  })
  # create list of lists w/ sorted gene names
  sorted.mat <- apply(pVal.mat, 1, function(x) {
    return(names(sort(x)))
  })
  # create candidate feature sets
  feat.sets <- list()
  for (i in 1:num.iter) {
    model.feats <- unique(unlist(as.list(sorted.mat[1:i,])))
    feat.sets[[i]] <- model.feats
  }
  return(feat.sets)
}

#' Performs feature selection using a pairwise wilcox test.
#' Each pairwise comparison is treated as a separate vector of features.
#' 
#' @param dat.mat Matrix of data (features X samples)
#' @param dat.labels Vector or factor of class labels
#' @param num.iter Number of candiate feature sets to return. Default of 25.
#' @return Lists of lists; each sub-list is a candidate feature set.
pwFeatsFull <- function(dat.mat, dat.labels, num.iter = 15) {
  dat.labels <- as.factor(dat.labels)
  dat.classes <- levels(dat.labels)
  # calculate p.values for each gene
  pVal.mat <- apply(dat.mat, 1, function(x) {
    p.mat <- pairwise.wilcox.test(x, dat.labels)$p.value
    p.list <- unlist(as.list(p.mat)); p.list <- p.list[!is.na(p.list)]
    return(p.list)
  })  
  # create matrix
  sorted.mat <- apply(pVal.mat, 1, function(x) {
    return(names(sort(x)))
  })
  # create candidate feature sets
  feat.sets <- list()
  for (i in 1:num.iter) {
    model.feats <- unique(unlist(as.list(sorted.mat[1:i,])))
    feat.sets[[i]] <- model.feats
  }
  return(feat.sets)
}

#' Performs feature selection using a pairwise wilcox test.
#' Each pairwise comparison is treated as a separate vector of features.
#' Travels down each list to get a new item from each comparison at each iteration.
#' 
#' @param dat.mat Matrix of data (features X samples)
#' @param dat.labels Vector or factor of class labels
#' @param num.iter Number of candiate feature sets to return. Default of 25.
#' @return Lists of lists; each sub-list is a candidate feature set.
pwFeatsFullUnique <- function(dat.mat, dat.labels, num.iter = 15) {
  dat.labels <- as.factor(dat.labels)
  dat.classes <- levels(dat.labels)
  # calculate p.values for each gene
  pVal.mat <- apply(dat.mat, 1, function(x) {
    p.mat <- pairwise.wilcox.test(x, dat.labels)$p.value
    p.list <- unlist(as.list(p.mat)); p.list <- p.list[!is.na(p.list)]
    return(p.list)
  })  
  # create matrix
  sorted.mat <- apply(pVal.mat, 1, function(x) {
    return(names(sort(x)))
  })
  # create candidate feature sets
  feat.sets <- list()
  for (i in 1:num.iter) {
    # create feature set 
    feat.set <- sorted.mat[1:i,1]
    for (j in 2:ncol(sorted.mat)) {
      iter.depth <- length(feat.set) + i
      new.feats <- setdiff(sorted.mat[1:iter.depth, j], feat.set)
      feat.set <- c(feat.set, new.feats[1:i])
    }
    # add to list
    feat.sets[[i]] <- feat.set
  }
  return(feat.sets)
}

#' Generates a train / test split of the given data
#' 
#' @param dat.mat Matrix of data (features X samples)
#' @param dat.labels Vector of class labels
#' @param train.per Percent of data to be used as training; default of 70%
#' @return A four items list; train.x (matrix), train.y (factor), test.x (matrix), test.y (factor)
TrainTestSplit <- function(dat.mat, dat.labels, train.per = 0.7) {
  # identify sample specific rates to perform an even split
  N <- ncol(dat.mat)
  sample.rates <- table(dat.labels); sample.rates <- sample.rates / sum(sample.rates)
  # select train / test samples
  train.samps <- c()
  for (cl in names(sample.rates)) {
    num.samps <- 0.7 * N * sample.rates[cl]
    train.samps <- c(train.samps, sample(which(dat.labels == cl), num.samps))
  }
  test.samps <- setdiff(1:N, train.samps)
  # generate train / test objects
  train.x <- dat.mat[, train.samps]; test.x <- dat.mat[, test.samps]
  train.y <- as.factor(dat.labels[train.samps]); test.y <- as.factor(dat.labels[test.samps])
  # return as list
  ret.list <- list('train.x' = train.x, 'train.y' = train.y, 'test.x' = test.x, 'test.y' = test.y)
  return(ret.list)
}

#' Performas multiclass AUROC analysis on a given model.
#' 
#' @param model Random Forest model object
#' @param test.x Matrix of test data (features X samples)
#' @param test.y Factor with class labels
#' @return Returns vector of aurocs for each class
mcAUROC <- function(model, test.x, test.y) {
  # make predictions
  test.predict <- predict(model, t(test.x), predict.all = TRUE)
  classes <- levels(test.y)
  ## get auroc for each class
  aurocs <- c()
  for (c.name in classes) {
    class.y <- as.character(test.y); class.y[which(class.y != c.name)] <- 'other'
    ntree <- ncol(test.predict$individual)
    class.percent <- apply(test.predict$individual, 1, function(x) {
      return(length(which(x == c.name)) / ntree)
    })
    class.roc <- roc(class.y, class.percent)
    aurocs <- c(aurocs, auc(class.roc))
  }
  ## return auroc for each class
  return(aurocs)
}

#' Generates a percent matrix from class votes of a random forest model
#'
#' @param predict.obj Prediction object.
#' @return Matrix of percentages for each class (samples X classes)
MakePercentMat <- function(predict.obj) {
  # get parameters
  ntree <- ncol(predict.obj$individual)
  class.labels <- levels(predict.obj$aggregate)
  samp.names <- rownames(predict.obj$individual)
  # get percentage vector
  percent.vect <- apply(predict.obj$individual, 1, function(x) {
    c.percent <- c()
    for (cl in class.labels) {
      c.percent <- c(c.percent, length(which(x == cl)) / ntree)
    }
    names(c.percent) <- class.labels
    return(c.percent)
  })
  # create matrix
  percent.mat <- matrix(percent.vect, nrow = 4)
  rownames(percent.mat) <- class.labels; colnames(percent.mat) <- samp.names
  return(t(percent.mat))
}

#' Returns the Youden's Index predictions for the ggiven data
#' 
#' @param percent.mat Percentage matrix from prediction.
#' @param yt.vect Named vector of Youden's Index threshold; names are class labels
#' @return Vector of YI classifications.
YIPredict <- function(percent.mat, yt.vect) {
  # get class labels
  class.labels <- names(yt.vect)
  # create yi matrix
  yi.mat <- percent.mat
  for (cl in class.labels) {
    yi.vect <- rep(0, nrow(yi.mat))
    yi.vect[which(yi.mat[,cl] >= yt.vect[cl])] <- 1
    yi.mat[,cl] <- yi.vect
  }
  # collapse yi matrix
  yi.class <- apply(yi.mat, 1, function(x) {
    if (sum(x) != 1) {
      return('unknown')
    } else {
      return(colnames(yi.mat)[which.max(x)])
    }
  })
  # return
  return(yi.class)
}

#' Generates and saves a conf table, normalized by colsums
#' 
#' @param pred.vect vector of predictions
#' @param class.vect vector of true classes
#' @param out.file name of file to write table to
WriteConfTable <- function(pred.vect, class.vect, out.file) {
  conf.table <- table(as.factor(pred.vect), as.factor(class.vect))
  conf.table <- round(t(t(conf.table) / colSums(conf.table)), digits = 3)
  write.csv(conf.table, file = out.file)
}

#' Generate color vector from percentage matrix.
#' 
#' @param percent.mat Percentage matrix from RF model.
#' @param class.colors Named list of colors for each class.
#' @param class.thresh Level at which to consider a class. Default of 0.2
#' @return Vector of mixed colors.
ClassThreshColVect <- function(percent.mat, class.colors, class.thresh = 0.2) {
  require(RColorBrewer)
  # generate color vector
  col.vect <- c()
  for (cell in rownames(percent.mat)) {
    # get highest two classes
    per.vect <- percent.mat[cell,]
    max.classes <- per.vect[names(sort(per.vect, decreasing = TRUE))[1:2]]
    # mix colors
    if (max.classes[2] > class.thresh) {
      max.classes <- max.classes / sum(max.classes)
      col.func <- colorRampPalette(class.colors[names(max.classes)])(11)
      col.val <- round(max.classes[1] * 10)
      col.vect <- c(col.vect, col.func[col.val])
    } else { # if second class doesn't pass threshold, don't mix
      col.vect <- c(col.vect, class.colors[names(max.classes)[1]])
    }
  }
  # name and return
  names(col.vect) <- rownames(percent.mat)
  return(col.vect)
}

#' Generates entropy vector from percentage matrix.
#' 
#' @param percent.mat Percentage matrix from RF model.
#' @return Vector of entropies for each sample
ClassEntropyVect <- function(percent.mat)  {
  entropy.vect <- c()
  for (s in rownames(percent.mat)) {
    samp.entropy <- -sum(sapply(percent.mat[s,], function(x) {
      if (x == 0) {return(0)}
      return(x * log(x, base = 2))
    }))
    entropy.vect <- c(entropy.vect, samp.entropy)
  }
  names(entropy.vect) <- rownames(percent.mat)
  return(entropy.vect)
}

#' Bootstrapped ISB log fold change analysis.
#' 
#' @param ref.mat Raw count matrix of reference samples
#' @param alt.mat Raw count matrix of alt / condition samples.
#' @param boot.num Number of botstraps. Default of 100.
#' @return Data frame w/ mean LFC between Alt and Ref, as well as SD of LFC in the bootstraps.
ISBLogFoldChange <- function(ref.mat, alt.mat, boot.num = 100) {
  # remove genes with zero counts
  ref.genes <- rownames(ref.mat)[which(rowSums(ref.mat) > 0)]
  alt.genes <- rownames(alt.mat)[which(rowSums(alt.mat) > 0)]
  shared.genes <- intersect(ref.genes, alt.genes)
  ref.mat <- ref.mat[shared.genes,]; alt.mat <- alt.mat[shared.genes,]
  # bootstrap LFC 
  lfc.mat <- matrix(0L, nrow = length(shared.genes), ncol = boot.num)
  rownames(lfc.mat) <- shared.genes
  for (b in 1:boot.num) {
    # generate bootsraps
    ref.boot <- rowSums(ref.mat[, sample(colnames(ref.mat), ncol(ref.mat), replace = TRUE)])
    alt.boot <- rowSums(alt.mat[, sample(colnames(alt.mat), ncol(alt.mat), replace = TRUE)])
    # normalize bootstraps
    ref.boot <- (ref.boot / sum(ref.boot)) * 1e6
    alt.boot <- (alt.boot / sum(alt.boot)) * 1e6
    # calculate lfc and store
    lfc.mat[,b] <- log(alt.boot / ref.boot, base = 2)
  }
  # compute mean / sd 
  mean.lfc <- rowMeans(lfc.mat); sd.lfc <- apply(lfc.mat, 1, sd)
  lfc.df <- data.frame('Mean' = mean.lfc, 'SD' = sd.lfc); rownames(lfc.df) <- shared.genes
  return(lfc.df)
}

gg_circle <- function(r, xc, yc, color="black", fill='NA', ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

#' Generates a circle plot for the given classification vector.
#' 
#' @param mvc.vect Vector of maximal vote classification.
#' @param percent.mat Matrix of vote percentages.
#' @param class.entropy Vector of clasification entropy.
#' @param plot.title TItle for this plot.
#' @param plot.path Path to save this plot.
#' @param class.cols Named vector of class colors.
#' @param p.vect Optional argument for protein activity vector. If included, used to color points.
#' @param p.name Optional argument for protein name. Used as title if included. Required if p.vect is provided.
#' @return Saves the plot object.
CirclePlot <- function(mvc.vect, percent.mat, class.entropy, plot.title, plot.path, class.cols, p.vect, p.name) {
  # plot parameters
  plot.rad <- 10; circle.int <- plot.rad / sqrt(2)
  text.buffer <- 1; anno.size <- 8
  max.ent <- log(4, base = 2)
  class.angles <- c(pi, 3*pi, 5*pi, 7*pi) / 4; names(class.angles) <- c('hsc', 'mpp', 'mlp', 'prog')
  percent.mat <- percent.mat[,names(class.angles)]
  # adjust entropy and calculate angles
  ent.vect <- (-1) * (class.entropy - max(class.entropy)) # invert so maximal entropy is at the origin
  ent.vect <- (plot.rad / max.ent) * ent.vect # rescale so that the maximum possible entropy would lie on the circle
  squared.per.mat <- percent.mat**2; squared.per.mat <- squared.per.mat / rowSums(squared.per.mat); samp.angles <- rowSums(t(t(squared.per.mat) * class.angles))
  # carteisan transformation
  x.coords <- ent.vect * cos(samp.angles)
  y.coords <- ent.vect * sin(samp.angles)
  if (!missing(p.vect)) {
    plot.dat <- data.frame('x.coord' = x.coords, 'y.coord' = y.coords)
    plot.dat[[p.name]] <- p.vect
    col.var <- p.name
  } else {
    plot.dat <- data.frame('x.coord' = x.coords, 'y.coord' = y.coords, 'class' = mvc.vect)
    col.var <- 'class'
  }
  # plot
  jpeg(file = plot.path,  width = 1100, height = 1050)
  circle.plot <- ggplot(plot.dat, aes(x.coord, y.coord)) + gg_circle(r = plot.rad, xc = 0, yc = 0, fill = 'white') + 
    geom_segment(x = -circle.int, y = -circle.int, xend = circle.int, yend = circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
    geom_segment(x = -circle.int, y = circle.int, xend = circle.int, yend = -circle.int, linetype = 'dashed', color = 'black', lwd = 1) + 
    xlim(c(-plot.rad, plot.rad)) + ylim(c(-plot.rad, plot.rad)) + ggtitle(plot.title) + 
    geom_point(aes_string(color = col.var)) + geom_density_2d(color = 'forestgreen') + 
    theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
          axis.title = element_blank(), axis.text = element_blank(), text = element_text(size = 16), panel.background = element_rect(fill = 'white', colour = 'white')) +
    annotate(geom="text", x = circle.int + text.buffer, y = circle.int + text.buffer, label="HSC", color = class.cols[['hsc']], fontface = 'bold', size = anno.size) +
    annotate(geom="text", x = -circle.int - text.buffer, y = circle.int + text.buffer, label="MPP", color = class.cols[['mpp']], fontface = 'bold', size = anno.size) +
    annotate(geom="text", x = -circle.int - text.buffer, y =  -circle.int - text.buffer, label="MLP", color = class.cols[['mlp']], fontface = 'bold', size = anno.size) +
    annotate(geom="text", x = circle.int + text.buffer, y = -circle.int - text.buffer, label="PROG", color = class.cols[['prog']], fontface = 'bold', size = anno.size)
  if (!missing(p.vect)) {
    circle.plot <- circle.plot + scale_color_gradient2(low = 'blue', high = 'red', mid = 'grey', midpoint = 0)
  } else {
    circle.plot <- circle.plot + scale_color_manual(values = class.cols)
  }
  print(circle.plot)
  dev.off()
}




