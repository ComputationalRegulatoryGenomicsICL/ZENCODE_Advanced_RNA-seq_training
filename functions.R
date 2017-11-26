

# In RNA-seq assay, parses the file with reads per gene

readTagsPerGene <- function(x){
  
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}

# Given the assay type, biosample and data pattern, lists all the files matching 
# pattern in the biosample

listFilesInBiosamples <- function(pathToData, assay, biosample, filesTypePattern){
  
  out <- list.files(path = paste(pathToData, assay, biosample, sep = "/"),
                    pattern = filesTypePattern, full.names = TRUE)
}

genericDESeq <- function(dataCounts, coldata, design, designVariations){
  head(dataCounts)
  dds <- DESeqDataSetFromMatrix(countData = dataCounts,
                                colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds$condition <- relevel(dds$condition, ref = designVariations[1])
  dds 
  dds <- DESeq(dds)
  res <- results(dds)
  summary(res)
  resLFC <- lfcShrink(dds, coef=2)
  summary(resLFC)
  plotMA(res, ylim =c(-2, 2))
  plotMA(resLFC, ylim=c(-2,2))
  d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                  returnData=TRUE)
  ggplot(d, aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(100,1000,5000))
  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  ntd <- normTransform(dds)
  meanSdPlot(assay(ntd))
  meanSdPlot(assay(vsd))
  meanSdPlot(assay(rld))
  plotPCA(vsd, intgroup=c("condition"))
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- vsd$condition
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  vsd
  
}

custom.see.genes <- function (data, edesign = data$edesign, time.col = 1, repl.col = 2, 
          group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
          cluster.data = 1, groups.vector = data$groups.vector, k = 9, 
          k.mclust = FALSE, cluster.method = "hclust", distance = "cor", 
          agglo.method = "ward.D", show.fit = FALSE, dis = NULL, step.method = "backward", 
          min.obs = 3, alfa = 0.05, nvar.correction = FALSE, show.lines = TRUE, 
          iter.max = 500, summary.mode = "median", color.mode = "rainbow", 
          cexlab = 1, legend = TRUE, newX11 = TRUE, ylim = NULL, main = NULL, 
          item = "genes", ...) 
{
  time = edesign[, time.col]
  repvect = edesign[, repl.col]
  groups = edesign[, group.cols]
  narrays <- length(time)
  if (!is.null(dim(data))) {
    dat <- as.data.frame(data)
    clusterdata <- data
  }
  else {
    clusterdata <- data[[cluster.data]]
    dat <- as.data.frame(data$sig.profiles)
  }
  if (nrow(dat) > 1) {
    dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
                                  1):ncol(dat)])
    count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
                       length(unique(repvect))), ]
    clusterdata <- dat
    if (any(is.na(clusterdata))) {
      if (cluster.method == "kmeans" || cluster.method == 
          "Mclust") {
        if (all(cluster.data != 1, cluster.data != "sig.profiles")) {
          clusterdata[is.na(clusterdata)] <- 0
        }
        else {
          mean.replic <- function(x) {
            tapply(as.numeric(x), repvect, mean, na.rm = TRUE)
          }
          MR <- t(apply(clusterdata, 1, mean.replic))
          if (any(is.na(MR))) {
            row.mean <- t(apply(MR, 1, mean, na.rm = TRUE))
            MRR <- matrix(row.mean, nrow(MR), ncol(MR))
            MR[is.na(MR)] <- MRR[is.na(MR)]
          }
          data.noNA <- matrix(NA, nrow(clusterdata), 
                              ncol(clusterdata))
          u.repvect <- unique(repvect)
          for (i in 1:nrow(clusterdata)) {
            for (j in 1:length(u.repvect)) {
              data.noNA[i, repvect == u.repvect[j]] = MR[i, 
                                                         u.repvect[j]]
            }
          }
          clusterdata <- data.noNA
        }
      }
    }
    if (!is.null(clusterdata)) {
      k <- min(k, nrow(dat), na.rm = TRUE)
      if (cluster.method == "hclust") {
        if (distance == "cor") {
          dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                            nrow(clusterdata), nrow(clusterdata)) - cor(t(clusterdata), 
                                                                        use = "pairwise.complete.obs")
          clust <- hclust(as.dist(dcorrel), method = agglo.method)
          c.algo.used = paste(cluster.method, "cor", 
                              agglo.method, sep = "_")
        }
        else {
          clust <- hclust(dist(clusterdata, method = distance), 
                          method = agglo.method)
          c.algo.used = paste(cluster.method, distance, 
                              agglo.method, sep = "_")
        }
        cut <- cutree(clust, k = k)
      }
      else if (cluster.method == "kmeans") {
        cut <- kmeans(clusterdata, k, iter.max)$cluster
        c.algo.used = paste("kmeans", k, iter.max, sep = "_")
      }
      else if (cluster.method == "Mclust") {
        if (k.mclust) {
          my.mclust <- Mclust(clusterdata)
          k = my.mclust$G
        }
        else {
          my.mclust <- Mclust(clusterdata, k)
        }
        cut <- my.mclust$class
        c.algo.used = paste("Mclust", k, sep = "_")
      }
      else stop("Invalid cluster algorithm")
      if (newX11) 
        X11()
      groups <- as.matrix(groups)
      colnames(groups) <- names.groups
      if (k <= 4) 
        par(mfrow = c(2, 2))
      else if (k <= 6) 
        par(mfrow = c(3, 2))
      else if (k > 6) 
        par(mfrow = c(3, 3))
      rn <- row.names(dat)
      dat1 <- dat
      for (x in rn){
        dat1[x, ] <- dat1[x, ] / max(dat1[x, ])
      }
      for (i in 1:(k)) {
        custom.PlotProfiles(data = dat1[cut == i, ], repvect = repvect, 
                     main = i, ylim = ylim, color.mode = color.mode, 
                     cond = rownames(edesign), item = item, ...)
      }
      if (newX11) 
        X11()
      if (k <= 4) {
        par(mfrow = c(2, 2))
        cexlab = 0.6
      }
      else if (k <= 6) {
        par(mfrow = c(3, 2))
        cexlab = 0.6
      }
      else if (k > 6) {
        par(mfrow = c(3, 3))
        cexlab = 0.35
      }
      for (j in 1:(k)) {
        PlotGroups(data = dat[cut == j, ], show.fit = show.fit, 
                   dis = dis, step.method = step.method, min.obs = min.obs, 
                   alfa = alfa, nvar.correction = nvar.correction, 
                   show.lines = show.lines, time = time, groups = groups, 
                   repvect = repvect, summary.mode = summary.mode, 
                   xlab = "time", main = paste("Cluster", j, sep = " "), 
                   ylim = ylim, cexlab = cexlab, legend = legend, 
                   groups.vector = groups.vector, item = item, 
                   ...)
      }
    }
    else {
      print("warning: impossible to compute hierarchical clustering")
      c.algo.used <- NULL
      cut <- 1
    }
  }
  else if (nrow(dat) == 1) {
    if (newX11) 
      X11()
    PlotProfiles(data = dat, repvect = repvect, main = NULL, 
                 ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
                 ...)
    if (newX11) 
      X11()
    PlotGroups(data = dat, show.fit = show.fit, dis = dis, 
               step.method = step.method, min.obs = min.obs, alfa = alfa, 
               nvar.correction = nvar.correction, show.lines = show.lines, 
               time = time, groups = groups, repvect = repvect, 
               summary.mode = summary.mode, xlab = "time", main = main, 
               ylim = ylim, cexlab = cexlab, legend = legend, groups.vector = groups.vector, 
               ...)
    c.algo.used <- NULL
    cut <- 1
  }
  else {
    print("warning: NULL data. No visualization possible")
    c.algo.used <- NULL
    cut <- NULL
  }
  OUTPUT <- list(cut, c.algo.used, groups)
  names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups")
  OUTPUT
}

custom.PlotProfiles <- function (data, cond, main = NULL, cex.xaxis = 0.5, ylim = NULL, 
          repvect, sub = NULL, color.mode = "rainbow", item = NULL) 
{
  pos.vline <- repvect[1:length(repvect)] - c(0, repvect[1:(length(repvect) - 
                                                              1)])
  for (i in 1:length(pos.vline)) {
    if (pos.vline[i] != 0) 
      pos.vline[i] <- i
    else pos.vline[i] <- pos.vline[i - 1]
  }
  if (is.null(ylim)) {
    ylim = c(min(as.matrix(data), na.rm = TRUE) * 1.1, max(as.matrix(data), 
                                                           na.rm = TRUE) * 1.1)
  }
  if (!is.vector(data)) {
    n = dim(data)[2]
    m = dim(data)[1]
    if (m == 1) 
      nom <- rownames(data)
    else nom <- NULL
    plot(x = c(1:n), y = data[1, ], type = "l", col = 1, 
         ylim = ylim, ylab = "expression value", xlab = " ", 
         main = paste("Cluster", main, "(", m, item, ")", 
                      nom, sep = " "), xaxt = "n")
    axis(1, at = 1:n, labels = substr(cond, 1, 26), cex.axis = cex.xaxis, 
         las = 2)
    if (color.mode == "rainbow") {
      abline(v = pos.vline, col = "light gray")
      for (i in 1:dim(data)[1]) {
        lines(x = c(1:n), y = data[i, ], col = i)
      }
    }
    else if (color.mode == "gray") {
      for (i in 1:dim(data)[1]) {
        lines(x = c(1:n), y = data[i, ], col = "black", alpha = 0.5)
      }
      yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
      lines(x = c(1:n), y = yy, col = "black")
    }
    else stop("Invalid mode, must be one of rainbow, gray")
  }
  else {
    n = length(data)
    plot(x = c(1:n), y = data, type = "l", col = 1, ylim = ylim, 
         ylab = "expression value", sub, xaxt = "n", xlab = " ")
    axis(1, at = 1:n, labels = cond, cex.axis = cex.xaxis, 
         las = 2)
    abline(v = pos.vline, col = "light gray")
  }
}