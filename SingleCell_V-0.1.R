## for single-cell sequencing transcriptome and 16sRNA sequencing
##
## 1. single-cell and 16s data usually has many '0's. These '0's do not mean it is not expressed. 
##    They do not even always mean very lowly expressed if the total read number (unique UMIs) is low. 
## 2. Especially for gut microbiota data (16s), as we all know there are 10^14 microbes in gut. so even 1e-5 
##    mean that there are 1 billion of the microbe. 
## 3. SOLUTION: 1) make '0's as NA. 2) use na.rm whenever fit. 3) use 'use="pair"' in cor function for clustering.
##

if(0){
  bioc.install("RcisTarget")
  bioc.install('GENIE3')
  bioc.install('AUCell')
  
  { ## SCENIC
    ## https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html
    install.packages("devtools")
    devtools::install_github("aertslab/SCopeLoomR")
    devtools::install_github("aertslab/SCENIC")
  }
  
  { ## SCUBA
    ## http://www.maths.bath.ac.uk/~masas/scuba/downloads.shtml
    ## R CMD INSTALL /hdd/software/SCUBA_Linux64_1.0.tar.gz
    
  }
  
  install.packages('Seurat')
}
source('~/codes/affymetrix/affy_survival_procedure.R')
## 1. QC
## 1.1. libsize.eff means library size effects.
##      d is a matrix, with columns as samples, rows as features.
zero2na <- function(d){
  d[d==0] <- NA
  invisible(d)
}

libsize.eff <- function(d){
  libsizes <- apply(d, 2, sum)
  detected.n <- apply(d, 2, function(x) sum(x>0))
  plot(libsizes, detected.n)
  my.plot.density(libsizes, main='libsize')
  my.plot.density(detected.n, main='number of detected features (>0)')

  d.na <- zero2na(d)  
  av.na <- apply(d.na, 1, mean, na.rm=T)
  ss.na <- apply(d.na, 1, sd, na.rm=T)
  plot(av.na, ss.na/av.na, main='after zero to NA')
  
  av <- apply(d, 1, mean)
  ss <- apply(d, 1, sd)
  plot(av, ss/av, main='before zero to NA')
  
  smoothScatter(log10(av), log10(av.na), main='mean')
  abline(0,1)
  smoothScatter(log10(ss), log10(ss.na), main='sd')
  abline(0,1)
  smoothScatter(ss/av, ss.na/av.na, main='sd/av')
  abline(0,1)
  ls.seg <- segmentize(d, libsizes, log='xy', main='libsize')
  dn.seg <- segmentize(d, detected.n, log='xy', main='detected')
  
  k <- ncol(ls.seg$with.na)
  smoothScatter(log10(ls.seg$with.na[,1]), log10(ls.seg$with.na[, k]), main="log10 libsize grouped")
  abline(0,1)
  smoothScatter(log10(ls.seg$with.0[,1]), log10(ls.seg$with.0[, k]), main='log10 libsize grouped')
  abline(0,1)

  smoothScatter(log10(dn.seg$with.na[,1]), log10(dn.seg$with.na[, k]), main="log10 detected.n grouped")
  abline(0,1)
  smoothScatter(log10(dn.seg$with.0[,1]), log10(dn.seg$with.0[, k]), main='log10 detected.n grouped')
  abline(0,1)
  
  invisible(list(libsize=libsizes, detected.n=detected.n, dn.seg=dn.seg, ls.seg=ls.seg, 
                 av=av, sd=ss, av.na=av.na, sd.na=ss.na))
}

mean.na <- function(x) mean(x, na.rm=TRUE)

## x is a type of statistics for each column of d, such as x = apply(d, 2, mean)
## d could have 0s, but not NAs.
segmentize <- function(d, x, compare.extreme=TRUE, pairs=FALSE, log=c('', 'xy')[1], main='', ...){
  probs <- c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1)
  x.cuts <- unique(quantile(x, probs=c(0, probs)))
  x.cuts[x.cuts == min(x.cuts)] <- min(x.cuts) - 1
  x.cuts[x.cuts == max(x.cuts)] <- max(x.cuts) + 1
  x.segs <- cut(x, x.cuts)
  names(x.segs) <- names(x)
  
  d.na <- d
  d.na[d.na == 0] <- NA
  d.segs.na <- aggregate.II(d.na, x.segs[colnames(d.na)], f='mean.na', byrow=FALSE)
  d.segs <- aggregate.II(d, x.segs[colnames(d)], f='mean', byrow=FALSE)
  
  if(compare.extreme){
    plot(d.segs.na[, 1], d.segs.na[, ncol(d.segs.na)], main=paste(main, ': with NA'), log=log); abline(0,1, col=2)
    plot(d.segs[, 1], d.segs[, ncol(d.segs)], main=paste(main, ': with 0'), log=log); abline(0,1, col=2)
  }
  if(pairs){
    pairs(d.segs, upper.panel=function(x, y){points(x, y);abline(0,1)}, main=paste(main, ': with 0'), log=log)
    pairs(d.segs.na, upper.panel=function(x, y){points(x, y);abline(0,1)}, main=paste(main, ': with NA'), log=log)
  }
  invisible(list(with.na=d.segs.na, with.0=d.segs))
}

sc.deg <- function(d, group1, group2, group1.name='group1', group2.name='group2', 
                   method = c('fisher', 'phyer', 'prop.test', "ttest", "wilcox", "sam")[1], 
                   log=FALSE, base=2, ...){
  if(method %in% c('fisher', 'phyper', 'binorm')){
    umi.size <- apply(d, 2, sum)
    group1.count <- apply(d[, group1], 1, sum)
    group2.count <- apply(d[, group2], 1, sum)
    raw.d <- cbind(group1.gc=group1.count, group2.gc=group2.count, 
                group1.tc = sum(umi.size[group1]), group2.tc = sum(umi.size[group2]),
                group1.size = length(group1), group2.size=length(group2))
    
    if(method == 'fisher'){
      p.l <- apply(raw.d, 1, function(x){
        wrapper.fisher.test(q=x[1], m=x[1] + x[2], n = x[3]+x[4]-x[1]-x[2], k = x[3],
                            alternative = c("two.sided"))
        })
      ps <- t(sapply(p.l, function(x) c(x$p.value, x$or)))
      colnames(ps) <- c('pvalue', 'or')
      ps <- ps[order(ps[, 'pvalue']), ]
      fdr <- p.adjust(ps[, 'pvalue'], method='fdr')
      res <- cbind(raw.d[rownames(ps), ], ps,  fdr=fdr)
    }else{
      if(method == 'phyper'){
        ps.t <- apply(raw.d, 1, function(x) phyper(x[1], x[1] + x[2], n = x[3]+x[4]-x[1]-x[2], k = x[3]))
      }
      if(method == 'prop.test'){
        ps.t <- apply(raw.d, 1, function(x) prop.test(x[1], x[3], p=x[2]/x[4])$p.value)  
      }
      ps <- 1 - abs((ps.t - 0.5))*2
      ps <- sort(ps)
      fdr <- p.adjust(ps, method='fdr')
      ## raw.d <- cbind(raw.d, apply(raw.d, 1, function(x) (x[1]/x[3])/(x[2]/x[4])))
      or <- (raw.d[,1]/raw.d[,3])/(raw.d[, 2]/raw.d[, 4])
      res <- cbind(raw.d[names(ps), ], pvalue=ps, or=or[names(ps)],  fdr=fdr)
    }
    res <- res[order(res[, 'pvalue']), ]
  }else{
    warning("d should have been normlized.")
    d[d == 0] <- NA
    if(log){
      d <- log(d, base=base) 
    }
    
    res <- deg(d, group1=group1, group2=group2, group1.name=group1.name, group2.name=group2.name, 
               method = method, log = FALSE, ...)  ## no log here for 'deg' function.
  }
  invisible(res)
}

sc.hclust <- function(d, log=FALSE, base=2, ...){
  d[d == 0] <- NA
  if(log){
    d <- log(d, base=base) 
  }

  res <- my.hclust(d, cor.use='pair',  ...)
  invisible(res)
}

sc.cor <- function(d, log=FALSE, base=2, ...){
  d[d==0] <- NA
  if(log){
    d <- log(d, base=base) 
  }

  res <- cor(d, use='pair', ...)
  invisible(res)
}

## For modeling purpose, use mean of non-0 to replace 0s or NAs. 
## remove rows, in whick a group of values are all 0s or NAs.
sc.impute <- function(d, sample.label, na.rm=TRUE){
  d[d==0] <- NA
  res <- NULL
  for(ss in split(names(sample.label), sample.label)){
    res <- cbind(res, t(apply(d[, ss], 1, function(x){
      if(sum(is.na(x)) < length(x)){
        x[is.na(x)] <- mean(x, na.rm=T)
        }
      x
      }
      )))
  }
  if(na.rm){
    res <- res[apply(res, 1, function(x) sum(is.na(x))) > 0, ]
  }
  invisible(res)
}

plot.mean.sd <- function(d, ...){
  av <- apply(d, 1, mean, na.rm=T)
  sd <- apply(d, 1, sd, na.rm=T)
  plot(av, sd, xlab='mean', ylab='sd', ...)
  plot(av, sd/av)
  invisible(sd/av)
}

## 'mean.n0' == 'mean without 0 (tranform 0 to NA)'.
## 'tmm', 'rle', 'upperquartile' come from edgeR. 
## TMM: weighted trimmed mean of M-values (to the reference.
## RLE: relative log expression. scale factor = median ratio of each sample to the median library.
## DESeq: scale factor = median ratio of each gene; ratio == x/geometric_mean(x) (x is gene signals across samples). 
## SCnorm is extremely slow. Be cautious!!!
sc.normalize <- function(d, method=c('libsize', 'mean.n0', 'quantile999', 'TMM', 'RLE', 
                                     'upperquartile', 'DESeq', 'scran', 'SCnorm')){
  method <- match.arg(method)
  
  nfs.norm <- function(d, nfs) t(t(d) * nfs)
  
  if(method=='libsize'){
    nfs <- 1e6/apply(d, 2, sum)
    d.norm <- nfs.norm(d, nfs)
  }

  if(method == 'mean.n0'){
    d[d==0] <- NA
    mean.n0 <- apply(d, 2, mean, na.rm=T)
    nfs <- mean(mean.n0)/mean.n0
    d.norm <- nfs.norm(d, nfs)
  }
  
  if(method == 'quantile999'){
    q999 <- apply(d, 2, quantile, probs = 0.999)
    nfs <- mean(q999)/q999
    d.norm <- nfs.norm(d, nfs)
  }

  if(is.element(method, c('TMM', 'tmm', 'rle', 'RLE', 'upperquartile'))){
    library(edgeR)
    nfs <- calcNormFactors(d, method=method)
    d.norm <- nfs.norm(d, nfs)
  }
  
  if(method == 'DESeq'){
    library(DESeq)
    nfs <- estimateSizeFactorsForMatrix(d)
    d.norm <- nfs.norm(d, nfs)
  }
  
  if(method=='scran'){
    library(scran)
    nfs <- 1/computeSumFactors(d)
    d.norm <- nfs.norm(d, nfs)
  }
  
  if(method=='SCnorm'){  ## extremely slow. Be cautious!!!
    library(SCnorm)
    DataNorm <- SCnorm(d, Conditions=rep(1, ncol(d)))
    d.norm <- DataNorm@metadata$NormalizedData
    nfs <- apply(d.norm, 2, sum)/apply(d, 2, sum)
  }
    
  invisible(list(method=method, nfs=nfs, d.norm=d.norm))
}

check.mito.gene.percentage <- function(d){
  mito.gi <- get.mito.genes()
  mito.geneids <- mito.gi[, "GeneID"]
  ## mito.geneids <- mito.gi[, "Symbol"]
  co.gs <- intersect(rownames(d), mito.geneids)
  percent.mito <- apply(d[co.gs, ], 2, sum)/apply(d, 2, sum)
  invisible(list(mito.genes = co.gs, percent.mito=percent.mito))
}

## source of cell cycle genes: https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8  Figure 4.
cell.cycle.genes <- c("6790", "9212", "891",  "9133", "55388",  "4171", "4172", "4173", "4174", 
                      "4175", "4176") 
names(cell.cycle.genes) <- c('AURKA', 'AURKB', 'CCNB1', 'CCNB2', 'MCM10', 'MCM2', 'MCM3', 
                             'MCM4', 'MCM5', 'MCM6', 'MCM7')
mito.genes <- get.mito.genes()

check.geneset <- function(d, norm.d, gs, gs.name, plot.each.gene=FALSE, cell.class.vector=NULL,
                          output.prefix, row.clust = F, column.clust = F, display.scale='row', 
                          h = NULL, cex=1, ...){
  cat("gs is genenames, names(gs) is geneID(entrezID)!\n")
  n0 <- length(gs)
  gs <- gs[is.element(gs, rownames(norm.d))]
  cat(n0, length(gs), "\n")
  
  gs.id2name <- names(gs)
  names(gs.id2name) <- gs
  
  if(max(norm.d, na.rm = T) > 30){
    norm.d <- log201(norm.d)
  }
  
  if(length(gs) < 3){
    warning('too few genes.')
    return()
  }
  
  gene.cc <- sc.cor(t(norm.d[gs,]), log=TRUE) 
  cell.cc <- NULL
  cell.cc <- sc.cor(d[gs,], log=TRUE)
  
  if(!is.null(cell.class.vector)){
    cell.class.vector <- sort(cell.class.vector)
    d <- d[, names(cell.class.vector)]
    norm.d <- norm.d[, names(cell.class.vector)]
  }else{
    cell.class.vector <- rep(1, ncol(d))
    names(cell.class.vector) <- colnames(d)
  }
  
  if(plot.each.gene){
    for(g in gs){
      ##my.plot.density(d[g,])
      barplot(table(norm.d[g, norm.d[g, ] > 0]), main=gs.id2name[g], log='y')
      my.plot.density(log2(norm.d[g, ]+1), main=gs.id2name[g])
    }
  }
  
  dh <- norm.d[gs, ]
  dh[is.na(dh)] <- 0
  ## dh <- log2(t(t(dh) + 1 + rnorm(ncol(norm.d), 0, 0.01)))
  rownames(dh) <- gs.id2name[rownames(dh)]
  ## l <- apply(dh, 2, function(x) length(unique(x)))
  ## dh <- dh[, l > 1]
  
  b <- apply(dh, 1, max) == 0
  cat("The following genes were not expressed in any cell: \n", rownames(dh)[b], "\n")
  dh <- dh[!b, ]
  
  if(is.null(h)){
    h <- 7 * (1 + nrow(dh)/100)    
  }
  
  if(nrow(dh) < 50){
    cex <- 2 * (1 + (50 - nrow(dh))/50)
  }
  pdf(file=paste(output.prefix, '_', gs.name, '.pdf', sep=''), width = 50, height = h)
  htmp.res <- my.heatmap(dh, column.clust = column.clust, row.clust = row.clust, 
                         col.center.zero = F, row.label = 'as.is', column.label = NULL,
                         col = color.gradient(low='white', high='red', n=100), X11 = F, 
                         column.class = split(names(cell.class.vector), cell.class.vector),
                         columnsep = T, sep.color = 'black', display.scale = display.scale, 
                         cex = cex, ...)
  dev.off()
  
  invisible(list(geneset=gs, gene.cc=gene.cc, cell.cc=cell.cc, heatmap.res=htmp.res, plot.height=h))
}  

check.single.gene <- function(d, norm.d, gene){
  if(!is.element(gene, rownames(d))){
    stop('gene is not in the experiment.')
  }

  x <- d[gene, ]
  print(summary(x))
  print(table(x > 0))
  barplot(c(expressed=sum(x>0), non.express=sum(x==0)), ylab='number of cells', main=gene)
  x1 <- x
  names(x1) <- NULL
  barplot(sort(x1), ylab='expression value', main=gene)
  umi.size <- apply(d, 2, sum)
  smoothScatter(log2(umi.size), log2(x+1), xlab='log2 total UMI', ylab='log2 gene expression', 
                main=gene)
  ## plot(density(log2(umi.size)[x > 0]), main=gene)
  ## lines(density(log2(umi.size)[x == 0]), col=2)
  my.vioplot(list(expressed=log2(umi.size)[x > 0], non.expressed=log2(umi.size)[x == 0]), 
             ylab='log2 total UMI', main=gene)
  cc <- sc.cor2(x, t(d))
  invisible(cc)
}

sc.cor2 <- function(x, d, log=FALSE, base=2, ...){
  n <- sum(x != 0 & !is.na(x))
  x[x==0] <- NA
  d[d==0] <- NA
  if(log){
    x <- log(x, base=base)
    d <- log(d, base=base)     
  }
  
  cc <- cor(x, d, use='pair', ...)
  
  used.ns <- apply(d, 2, function(x, y) sum(!is.na(x) & !is.na(y)), y=x)
  cc <- cbind(cc=cc[1, ], used.n=used.ns)
  cc <- cc[order(cc[,1], decreasing = TRUE), ]
  invisible(list(n=n, total.n=length(x), cc=cc, b5=cc[, 'used.n'] >= 5))
}
## For a pair of Genes (g1 and g2), who are correlated at both expressed or some specifical points, 
## but do not know whether others pairs (with only 1 expressed data point) are correlated or not.
## This function can only help to identfiy non-correleated pairs, but not potential correlated pairs,
## due to missing data.
## This function also not working on the pairs with no one point expressed.
identify.non.cor <- function(d, norm.d, g1, g2, norm.d.log=TRUE){
  ## identify.non.cor <- function(x, y, x0, y0){
  if(all.equal(rownames(d), rownames(norm.d)) & all.equal(colnames(d), colnames(norm.d))){
    ## do nothing.
  }else{
    stop('rownames or colnames of d and norm.d are not exactly same.')
  }
  
  if(norm.d.log){
    x <- log2(norm.d[g1, ] + 1)
    y <- log2(norm.d[g2, ] + 1)
  }else{
    x <- norm.d[g1, ]
    y <- norm.d[g2, ]
  }
  
  x0 <- d[g1,]
  y0 <- d[g2, ]
  
  plot(x, y)
  x.expressed <- x > 0 & !is.na(x)
  y.expressed <- y > 0 & !is.na(y)
  both.expressed <- x.expressed & y.expressed
  both.nonexpressed <- (!x.expressed) & (!y.expressed)
  x.expressed.only <- x.expressed & (!y.expressed)
  y.expressed.only <- y.expressed & (!x.expressed)
  
  if(sum(both.expressed) < 10){
    stop('less than 10 useful data pairs.')
  }
  
  t1 <- loess(x[both.expressed] ~ y[both.expressed])
  x[y.expressed.only] <- predict(t1, y[y.expressed.only])
  
  t2 <- loess(y[both.expressed] ~ x[both.expressed])
  y[x.expressed.only] <- predict(t2, x[x.expressed.only])
  lo.x <- x[both.expressed]
  lo.y <- t2$fitted
  o <- order(lo.x)
  lines(lo.x[o], lo.y[o], col=2, lty=2)
  b <- x.expressed.only
  points(x[b], y[b], pch='*', col=2)
  b <- y.expressed.only
  points(x[b], y[b], pch='*', col=3)
  
  t.f <- function(v, v0, x){
    xv <- cbind(x, v)
    t <- sort(xv)[which(order(xv) == 1) - 1]
    res <- v0[v == t][1]
    res
  }
  {
    m <- cbind(oo = 1:length(x), x=x, y=y, x0=x0, y0=y0, x1=x0, y1=y0, p=NA) ## oo == origin order
    ## m <- m[order(x), ]
    for(i in 1:nrow(m)){
      print(i)
      s <- rownames(m)[i]
      if(m[i, 'x'] != 0 & m[i, 'x0'] == 0 & !is.na(m[i, 'x'])){
        m[i, 'x1'] <- t.f(norm.d[, s], d[, s], m[i, 'x'])
        m[i, 'p'] <- pbinom(0, sum(x0), m[i, 'x1']/sum(x0))
      }
      
      if(m[i, 'y'] != 0 & m[i, 'y0'] == 0 & !is.na(m[i, 'y'])){
        m[i, 'y1'] <- t.f(norm.d[, s], d[, s], m[i, 'y'])
        m[i, 'p'] <- pbinom(0, sum(y0), m[i, 'y1']/sum(y0))
      }
    }
  }
  
  ## Small p means that at that point, the g1 and g2 did not follow correlation pattern shown in 
  ## both expressed cases.
  invisible(list(x=x, y=y, m=m))
}

sc.qc <- function(d, outFile){
  q <- apply(d, 2, quantile, 0.99)
  s <- apply(d, 2, sum)
  m <- apply(d, 2, max)
  
  pdf(outFile)
  
  my.plot.density(q, main='quantile0.99')
  my.plot.density(s, main='sum')
  my.plot.density(m, main='max')
  
  smoothScatter(q,s, xlab='quantile0.99', ylab='sum')
  smoothScatter(q,m, xlab='quantile0.99', ylab='max')
  smoothScatter(s,m, xlab='sum', ylab='max')
  
  tn <- apply(d, 2, sum)
  gn <- apply(d > 0, 2, sum)
  smoothScatter(tn, gn, xlab='total UMI', ylab='number of detected genes (reads > 0)', 
                ylim=c(0, max(gn)), xlim=c(0, max(tn)))
  abline(h=min(gn), v=min(tn), lty=2, col='grey')
  my.plot.density(gn, main='number of detected genes (>0)')
  
  gn.median <- median(gn)
  gn.mad <- mad(gn)
  gn.cut <- (gn.median - 3.1 * gn.mad)  ## qnorm(0.999) == 3.1s
  abline(v=gn.cut)
  
  mito <- check.mito.gene.percentage(d)
  my.plot.density(mito$percent.mito)
  mito.median <- median(mito$percent.mito)
  mito.mad <- mad(mito$percent.mito)
  mito.cut <- (mito.median + 3.1 * mito.mad) ## qnorm(0.001) == -3.1
  abline(v=mito.cut)
  print(table(mito$percent.mito > mito.cut))
  
  dev.off()
  
  invisible(list(mito=mito, mito.cut=mito.cut, good.size=tn, detected.gene.n=gn, detected.gene.cut = gn.cut))
}

## hard core analysis
## d is raw data with count
## norm.d is normalized d, no matter which way to be normalized
## cell.class.vector supposed to be 
## sampleName is the name of the singleCell Seq sample.
sc.hc <- function(d, norm.d, cell.class.vector, sampleName='this', gsea.b=TRUE, print.gsea=FALSE,
                  gsoi.l){
  if(!setequal(colnames(d), names(cell.class.vector))){
    stop(" colnames(d) != names(cell.class.vector). ")
  }

  load('~/data/geneAnnotation/human/hs.gi.RData')
  source('~/codes/pathway/my_gsea.R')

  ccv <- cell.class.vector
  ##if(sum(ccv[-1] != ccv[-length(ccv)]) > length(unique(ccv))){
  ccv <- sort(ccv)
  ## }

  ## tv <- table(ccv)
  ## ccv <- ccv[is.element(ccv, names(tv)[tv < 3])]
  
  d <- d[, names(ccv)]
  norm.d <- norm.d[, names(ccv)]
  cell.class.l <- split(names(ccv), ccv)
  
  d.grps <- aggregate.II(d, ccv[colnames(d)], f = 'mean', byrow = F)
  
  nc <- length(cell.class.l)
  load('~/data/geneAnnotation/human/ncbi.gene.info.RData')
  
  deg.l <- NULL
  for(i in 1:nc){
    if(length(cell.class.l[[i]]) < 3){
      next
    }
    deg.l[[as.character(i)]] <- sc.deg(d, group1 = cell.class.l[[i]], 
                                       group2 = unlist(cell.class.l[-i]), 
                                       group1.name = i, group2.name = 'other',
                                       method = 'phyper')
  }
  
  gsea.kegg.l <- gsea.go.l <- NULL
  
  if(gsea.b){
    for(i in 1:length(deg.l)){
      gsea.kegg <- batch.gsea(sort(deg.l[[i]][, 'or']), genesets = hs.gi$kegg)
      gsea.kegg <- data.frame(gsea.kegg, hs.gi$kegg.info[rownames(gsea.kegg), ])
      b <- gsea.kegg[, "pvalue"] < 0.05/nrow(gsea.kegg) & gsea.kegg[, "ES"] > 0
      gsea.kegg.l[[i]] <- gsea.kegg
      
      gsea.go <- batch.gsea(sort(deg.l[[i]][, 'or']), genesets = hs.gi$go)
      gsea.go <- data.frame(gsea.go, hs.gi$go.info[rownames(gsea.go), ])
      b <- gsea.go[, "pvalue"] < 0.05/nrow(gsea.go) & gsea.go[, "ES"] > 0 & gsea.go[, "ontology"] == 'BP'
      gsea.go.l[[i]] <- gsea.go
      
      if(print.gsea){
        print.n <- 30
        print(hs.gi$kegg.info[rownames(gsea.kegg)[b], , drop=F][1:min(sum(b), print.n), ])
        print(hs.gi$go.info[rownames(gsea.go)[b], 1, drop=F][1:min(sum(b), print.n), ])
      }
    }
    save(deg.l, gsea.go.l, gsea.kegg.l, file='middle.res.RData')
  }
  
  ## load('middle.res.RData')
  
  if(max(norm.d, na.rm = T) > 30){
    norm.d <- log201(norm.d)
  }
  
  marker.gene.l <- NULL
  ## my.hclust(raw.d[1:1000, ], sample.class = split(names(hc.10), hc.10))
  for(i in 1:length(deg.l)){
    ts <- deg.l[[i]]
    b.up <- d.grps[rownames(ts), as.character(i)] == apply(d.grps[rownames(ts), ], 1, max)
    b.down <- d.grps[rownames(ts), as.character(i)] == apply(d.grps[rownames(ts), ], 1, min)
    
    ts.up <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (ts[, 'or'] > 2 | is.infinite(ts[,'or'])) & b.up, ]
    ts.down <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (ts[, 'or'] < 0.5 & !is.infinite(ts[,'or'])) & b.down, ]
    
    co.up.gs <- intersect(rownames(gene.info$gene.info), row.names(ts.up))
    co.down.gs <- intersect(rownames(gene.info$gene.info), row.names(ts.down))
    co.gs <- union(co.up.gs, co.down.gs)
    marker.gene.l[[i]] <- data.frame(gene.info$gene.info[co.gs, ], ts[co.gs, ])
    
    td <- norm.d[co.gs, ]
    b <- ccv[colnames(td)] == i
    colnames(td)[b] <- '1'
    colnames(td)[!b] <- '.'
    td[is.na(td)] <- 0

    if(nrow(td) > 200){
      td <- td[1:200, ]
    }
    
    width=1200
    height = max(1000 * nrow(td)/2000, 500)
    tiff(file=paste(sampleName, '_heatmap-', i, '.tif', sep = ''), width = width, height = height)
    my.heatmap(td, row.clust = F, column.clust = F, 
               col.center.zero = F, column.label = 'as.is', row.label = NULL,
               col = color.gradient(low='white', high='red', n=100), 
               columnsep = T, sep.color = 'black',
               column.class = cell.class.l, X11 = F, main=i)
    dev.off()
  }
  
  all.markers <- sort(unique(unlist(sapply(marker.gene.l, rownames))))
  tiff(file=paste(sampleName, '_heatmap-allMarker.tif', sep=''), height=1000, width = 1200)
  td <- norm.d[all.markers, ]
  td[is.na(td)] <- 0
  my.heatmap(td, row.clust = T, column.clust = F, 
             col.center.zero = F, column.label = NULL, 
             col = color.gradient(low='white', high='red', n=100), 
             column.class = split(names(ccv), ccv), X11 = F)
  dev.off()
  
  invisible(list(deg.l=deg.l, marker.gene.l=marker.gene.l, gsea.kegg.l=gsea.kegg.l, gsea.go.l=gsea.go.l))
}

## 
sc.classification <- function(d, method=c('hclust', 'kmeans'), k){
  if(method == 'hclust'){
    origin.res <- sc.hclust(raw.d, noPlot=TRUE, noX = TRUE)
    class.v <- cutree(origin.res, k=k)
  }
  if(method == 'kmeans'){
    origin.res <- kmeans(t(raw.d), centers = 10)
    class.v <- origin.res$cluster
  }
  
  invisible(list(class.v=class.v, origin.res=origin.res, method=method, k=k))
}

check.geneset.batch <- function(d, norm.d, gsoi.l, min.sample.n=3, class.v, output.prefix='this', ...){
  res.l <- list()
  for(gs.name in names(gsoi.l)){
    gs <- gsoi.l[[gs.name]]
    x1 <- raw.d[intersect(rownames(raw.d), names(gs)), , drop=F]
    b <- apply(x1 > 0, 1, sum) > min.sample.n
    raw.d.t <- x1[b, ,drop=F]
    norm.d.t <- norm.d[rownames(raw.d.t), , drop=F]
    res.l[[gs.name]] <- check.geneset(raw.d.t, norm.d.t, gs, gs.name = gs.name, class.v = class.v, 
                                      output.prefix = output.prefix, ...)
  }
  invisible(res.l)
}

## data dir is 10x's data directory containing: barcodes.tsv, genes.tsv, matrix.mtx
read.10x.mtx <- function(data.dir, min.genes = 1000, ...){
  library(Seurat)
  t <- Read10X(data.dir = data.dir)
  t2 <- CreateSeuratObject(t, min.genes = min.genes, ...)
  t3 <- as.matrix(t2@raw.data)
  invisible(t3)
}

## d is output of read.10x.mtx
sc.preprocess <- function(d, min.expr.genes = 1000){
  dnames <- get.gene.info(rownames(d))
  d.entrez <- dm.rowname.convert(d, id.m = dnames[, c("alias", "GeneID")], keep.dup = TRUE)
  
  rm(d)
  mito=check.mito.gene.percentage(d.entrez)
  mito.median <- median(mito$percent.mito)
  mito.mad <- mad(mito$percent.mito)
  mito.cut <- (mito.median + 2* mito.mad)
  
  rd <- d.entrez[, names(mito$percent.mito)[mito$percent.mito < mito.cut]]
  rd <- rd[, apply(rd > 0, 2, sum) > min.expr.genes]
  norm.d <- sc.normalize(rd, method = 'mean.n0')
  
  mito.para = list(mito=mito, mito.cut=mito.cut)
  invisible(list(raw.d=rd, norm.d=norm.d, mito.para=mito.para, min.expr.genes=min.expr.genes))
}

## a QC test for classification/clustering results.
test.classes.libsize <- function(raw.d, class.v){
  this.d <- cbind(lib.size=apply(raw.d, 2, sum), class=class.v[colnames(raw.d)])
  res <- anova.simple.wrapper(this.d, y.column = 'lib.size', x.column = 'class', log = F, verbose = T)
  my.vioplot(split(this.d[, 1], this.d[,2]))
  invisible(res)
}

## define gene sets of interest
compile.geneset.of.interest <- function(){
  load('~/data/geneAnnotation/human/gencode.v27.annotation.RData')
  trav.ind <- grep('^TRAV', gencode$gene.gr$gene_name)
  t1 <- gencode$gene.gr$gene_name[trav.ind]
  t2 <- get.gene.info(t1)
  tt <- t2[, c('Symbol', 'GeneID')]
  tt <- unique.matrix(tt)
  trav.genes <- tt[, 'Symbol']
  names(trav.genes) <- tt[, 'GeneID']
  
  ## Human endometrial stromal cells: CD146，IGFBP-1, PRL, TF, PAI-1，CD45
  t <- c('CD146', 'IGFBP1', 'PRL', 'PAI-1', 'CD45') ## 'TF' not sure which genes
  tt <- get.gene.info(t)
  stromal.genes <- tt[, 'Symbol']
  names(stromal.genes) <- tt[, 'GeneID']
  
  ## Human endometrial epithelial cells:EpCAM
  t <- 'EPCAM'
  tt <- get.gene.info(t)
  epithelial.genes <- tt[, 'Symbol']
  names(epithelial.genes) <- tt[, 'GeneID']
  
  ## immuno genes
  t <- read.table('~/XKYBio/A_研发/Platforms/知识库/CellTypeMarkers/ImmuneAssociatedGenes.txt', 
                  sep='\t', header=T, as.is =T)
  t <- unique.matrix(t[, 2:3])
  immune.genes <- t[, 1]
  names(immune.genes) <- t[, 2]
  
  ## TFs
  load('~/data/geneAnnotation/human/TFs/human_tfs.RData')
  tfs <- tfs.info[, "Symbol"]
  names(tfs) <- tfs.info[, "GeneID"]
  
  gsoi.l <- list(immune=immune.genes, epithelial=epithelial.genes, stromal=stromal.genes, trav=trav.genes, 
                 tfs = tfs)
  
  immune.core <- gsoi.l[['immune']][is.element(gsoi.l[['immune']], 
                                               c('CD4', 'CD8A', 'CD8B', 'CD3G', 'CD3E', 'CD3D', 'CD19', 'PTPRC', 
                                                 'MS4A1', 'FOXP3', 'IFNG', 'PRF1', 'GZMA', 'GZMK', 'CTLA4', 'IDO1', 
                                                 'LAG3', 'TIM3', 'PDCD1', 'CD274', 'CD276', 'CCR7', 'HLA-A', 'HLA-B', 
                                                 'HLA-C', 'HLA-DRB1', 'HLA-DPB1', 'HLA-DQB1', 'TAP1', 'TAP2', 'B2M'))]

  gsoi.l[['immune.core']] <- immune.core

  t1 <- c('CD247', 'CD3G', 'CD3E', 'CD3D', ## CD3 family, cD247 = CD3-ZETA/CD3H/CD3Q/CD3Z 
         'FOXP3', ## Treg
         'CD4', 'CD8A', 'CD8B', 
         'CD34', ## stem cell, endothelial cell
         'CD19', 'MS4A1',   ## MS4A1 == CD20. B cell 
         'ITGAX', 'IL3RA',  ## ITGAX == CD11c; CD123==IL3RA. DC cell
         'NCAM1',   ## CD56 == NCAM1. NK cell
          'CD14', 'CD33', ## Macrophage/Monocyte
         'CEACAM8',   ## CD66b == CEACAM8. Granulocyte
         'ITGA2B', 'ITGB3', 'SELP',   ## CD41==ITGA2B, CD61== ITGB3, CD62 = SELP. Platelet
         'GYPA',   ## CD235a == GYPA. Erythrocyte
         'MCAM',  ## MCAM == CD146. Endothelial Cell
         'EPCAM'  ## CD326 == EPCAM. Epithelial Cell
         )
  t2 <- get.gene.info(t1)
  gsoi.l[['immune.cell.marker']] <- get.my.geneset.format(t2[,3])
  for(gs.name in names(gsoi.l)){
    gsoi.l[[gs.name]] <- get.my.geneset.format(gsoi.l[[gs.name]])
  }
  save(gsoi.l, file='~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
  invisible(gsoi.l)
}
