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
    ## For human
    dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
                 "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
    for(featherURL in dbFiles)
    {
      download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
      descrURL <- gsub(".feather$", ".descr", featherURL)
      if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
    }
    
    # For mouse
    ## dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
    ##              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
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
        ## ps.t <- apply(raw.d, 1, function(x) phyper(x[1], x[1] + x[2], n = x[3]+x[4]-x[1]-x[2], k = min(x[3:4])))
        if(raw.d[1,3] > raw.d[1,4]){
          ps.t <- apply(raw.d, 1, function(x) phyper(x[2], x[1], n = x[3] - x[1], k = x[4]))
        }else{
          ps.t <- apply(raw.d, 1, function(x) phyper(x[1], x[2], n = x[4] - x[2], k = x[3]))
        }
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
    b <- res[, 1] == 0 & res[,2] == 0
    res[b, 'pvalue'] <- NA
    res[b, 'fdr'] <- NA
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

check.single.gene.01 <- function(umi.d, gene=NULL, cell.sizes, simple.res=FALSE){
  if(is.matrix(umi.d)){
    ss <- colnames(umi.d)
    gi <- get.gene.info(gene)
    if(nrow(gi) > 1){
      print(gi[, 1:4])
      stop('please select one gene (GeneID can be unique).')
    }
    count <- umi.d[gi[1, 'GeneID'], ]
  }else{
    ss <- names(umi.d)
    count <- umi.d
  }
  
  if(length(ss) == 0 || length(setdiff(ss, names(cell.sizes))) > 0){
    stop('umi.d does not have name or cell.sizes is not complete.')
  }
  
  cell.sizes <- cell.sizes[names(count)]

  r <- count/cell.sizes
  min.r <- min(r[r > 0])
  test.res <- apply(cbind(count, cell.sizes), 1, 
                    function(x, p) prop.test(x[1], x[2], p=p, alternative = 'less'), 
                    p=min.r)
  ps <- sapply(test.res, function(x) x$p.value)
  res <- cbind(count=count, cell.size=cell.sizes, pvalue=ps, min.r=min.r)
  
  if(simple.res){
    invisible(ps)
  }else{
    invisible(res)
  }
}

get.present.m <- function(umi.d){
  res.m <- umi.d - umi.d
  cell.sizes <- apply(umi.d, 2, sum)
  d <- umi.d[apply(umi.d, 1, max) > 0, ]
  pm <- apply(d, 1, check.single.gene.01, cell.sizes = cell.sizes, simple.res=TRUE)
  pm <- t(pm)
  res.m[rownames(pm), colnames(pm)] <- pm
  res.m <- as.matrix(res.m)
  invisible(res.m)
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
  if(length(co.gs) == 0){
    percent.mito <- NULL
  }else{
    percent.mito <- apply(d[co.gs, ,drop=F], 2, sum)/apply(d, 2, sum)    
  }

  invisible(list(mito.genes = co.gs, percent.mito=percent.mito))
}

## source of cell cycle genes: https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8  Figure 4.
cell.cycle.genes <- c("6790", "9212", "891",  "9133", "55388",  "4171", "4172", "4173", "4174", 
                      "4175", "4176", '4288') 
names(cell.cycle.genes) <- c('AURKA', 'AURKB', 'CCNB1', 'CCNB2', 'MCM10', 'MCM2', 'MCM3', 
                             'MCM4', 'MCM5', 'MCM6', 'MCM7', 'MKI67')
mito.genes <- get.mito.genes()

check.geneset <- function(d, norm.d, gs, gs.name, plot.each.gene=FALSE, cell.class.vector=NULL,
                          output.prefix=NULL, row.clust = F, column.clust = F, display.scale='none', 
                          h = NULL, cex=1, main=NULL, x11.bulk=TRUE, ...){
  
  cat("gs is genenames, names(gs) is geneID(entrezID)!\n")
  n0 <- length(gs)
  gs <- gs[is.element(gs, rownames(norm.d))]
  cat(n0, length(gs), "\n")
  
  gs.id2name <- names(gs)
  names(gs.id2name) <- gs
  
  d.aggr <- aggregate.II(d, v=cell.class.vector, f = 'sum', byrow = F)
  d.aggr <- apply(d.aggr, 2, function(x) x/sum(x)*1e6)
  
  d <- d[gs, ]
  norm.d <- norm.d[gs, ]
  
  if(max(norm.d, na.rm = T) > 30){
    norm.d <- log201(norm.d)
  }
  
  if(length(gs) < 3){
    warning('too few genes.')
    return()
  }
  
  gene.cc <- sc.cor(t(norm.d[gs,]), log=FALSE) 
  cell.cc <- NULL
  if(column.clust){
    cell.cc <- sc.cor(d[gs,], log=FALSE)    
  }

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
      barplot(table(norm.d[g, norm.d[g, ] > 0]), main=paste(main,gs.id2name[g], sep='::'), log='y')
      my.plot.density(log2(norm.d[g, ]+1), main=paste(main, gs.id2name[g], sep='::'))
    }
  }
  
  dh <- norm.d[gs, ,drop=F]
  dh[is.na(dh)] <- 0
  ## dh <- log2(t(t(dh) + 1 + rnorm(ncol(norm.d), 0, 0.01)))
  rownames(dh) <- gs.id2name[rownames(dh)]
  ## l <- apply(dh, 2, function(x) length(unique(x)))
  ## dh <- dh[, l > 1]
  
  b <- apply(dh, 1, max) == 0
  cat("The following genes were not expressed in any cell: \n", rownames(dh)[b], "\n")
  dh <- dh[!b, ,drop=F]
  
  if(is.null(h)){
    h <- 7 * (1 + nrow(dh)/100)    
  }
  
  if(nrow(dh) < 50 & cex == 1){
    cex <- 2 * (1 + (50 - nrow(dh))/50)
  }
  
  main <- paste(main, gs.name, sep='::')
  
  ## d.aggr <- aggregate.II(add.noise(dh), v=cell.class.vector, f = 'mean', byrow = F)
  rownames(d.aggr) <- gs.id2name[rownames(d.aggr)]
  d.aggr <- d.aggr[rownames(dh), ,drop=F]
  
  if(x11.bulk){ ## plot in X11, and will be plot again in PDF.
    my.heatmap(d.aggr, column.clust = T, row.clust = T, col.center.zero = F, row.label = 'as.is', 
               column.label = 'as.is', col = color.gradient(low='white', middle = 'blue', high='red', n=100),
               X11 = T,
               display.scale = display.scale, grid = T, grid.color = 'grey', cex=2, 
               main=main)
  }
  
  htmp.res <- NA
  if(!is.null(output.prefix)){
    pdf(file=paste(output.prefix, '_', gs.name, '.pdf', sep=''), width = 50, height = h)
    on.exit(dev.off())
    x11.b <- FALSE
  }else{
    x11(width = 50, height = h)
    x11.b <- TRUE
  }
  
  if(nrow(dh) > 0){
    if(nrow(dh) == 1){
      row.clust <- FALSE
    }
    htmp.res <- my.heatmap(dh, column.clust = F, row.clust = row.clust, 
                           col.center.zero = F, row.label = 'as.is', column.label = NULL,
                           col = color.gradient(low='white', middle = 'blue', high='red', n=100), X11 = x11.b, 
                           column.class = split(names(cell.class.vector), cell.class.vector),
                           columnsep = T, sep.color = 'black', display.scale = display.scale, 
                           cex = cex, main=main, ...)
    if(column.clust & nrow(dh) > 1){  ## is column.clust is TRUE, still draw a no column clust plot}
      htmp.res <- my.heatmap(add.noise(dh, 0, 1e-10), column.clust = column.clust, row.clust = row.clust,   
                             col.center.zero = F, row.label = 'as.is', column.label = NULL,
                             col = color.gradient(low='white', middle = 'blue', high='red', n=100), X11 = x11.b, 
                             column.class = split(names(cell.class.vector), cell.class.vector),
                             columnsep = T, sep.color = 'black', display.scale = display.scale, 
                             cex = cex, main=main, ...)  
    }
  }

  { ## plot again in PDF
    my.heatmap(d.aggr, column.clust = T, row.clust = T, col.center.zero = F, row.label = 'as.is', 
               column.label = 'as.is', col = color.gradient(low='white', middle = 'blue', high='red', n=100), X11 = x11.b,
               display.scale = display.scale, grid = T, grid.color = 'grey', cex=2,
               main=main)
    
  }

  invisible(list(geneset=gs, gene.cc=gene.cc, cell.cc=cell.cc, heatmap.res=htmp.res, d.aggr=d.aggr))
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

plot.tsne <- function(tsne.d, cls, fast=TRUE, ...){
  ## tsne.file <- file.path(res_dir, sampleName, 'analysis/tsne/2_components/projection.csv')
  ## tsne <- read.table(tsne.file, as.is=T, header=T, row.names = 1, sep=',')
  ## plot(tsne[,1:2], pch='.')
  ## smoothScatter(t[,1:2], pch='*')
  
  ## cls.file <- file.path(res_dir, sampleName, 'analysis/clustering/graphclust/clusters.csv')
  ## cls <- read.table(cls.file, as.is=T, header=T, row.names = 1, sep=',')
  plot(tsne.d[,1:2], pch='.', ...)
  cols <- my.colors

  for(i in 1:max(cls)){
    if(fast){
      points(tsne.d[names(cls)[cls == i], 1:2], col=cols[(i-1)%%10+1], pch='*')
    }else{
      points(tsne.d[names(cls)[cls == i], 1:2], col=cols[(i-1)%%10+1], pch=(i-1)%%8+1)  
    }
  }
}

plot.tsne.II <- function(tsne.d, class.l, colStyle=c('rainbow', 'heat')[1], pch='.', 
                      cex=4, plot.legend=TRUE, main=NULL, ...){
  l <- length(class.l)
  plot(tsne.d, pch='.', main=main, ...)
  if(colStyle=='rainbow'){
    cols <- rainbow(l)      
  }
  if(colStyle == 'heat'){
    cols <- heat.colors(l)
  }
  if(l == 2){
    cols <- c('grey', 'red')
  }
  
  for(i in 1:l){
    points(tsne.d[class.l[[i]], ], col=cols[i], pch=pch, cex=cex)
  }
  if(plot.legend){
    if(l > 10){
      x11(height = l + 1) 
      plot(0,0, col='white', main=main)
      print('Please click a position to draw legend.')
      locator(1)
      legend('topright', col=cols, legend = names(class.l), pch='.', pt.cex=15, cex=1.5)
    }else{
      plot(0,0, col='white', main=main)
      legend('topright', col=cols, legend = names(class.l), pch='.', pt.cex=10, cex=1)
    }
  }
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
                  gsoi.l, output.dir='.'){
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
  
  if(max(norm.d, na.rm = T) > 30){
    norm.d <- log201(norm.d)
  }
  
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
                                       method = 'ttest')
  }
  
  gsea.kegg.l <- gsea.go.l <- NULL
  
  if(gsea.b){
    for(i in 1:length(deg.l)){
      gsea.res <- batch.gsea.count.data(count.data.deg = deg.l[[i]])
      gsea.kegg.l[[i]] <- gsea.res$kegg
      gsea.go.l[[i]] <- gsea.res[c('go.bp', 'go.mf', 'go.cc')]
      
      if(0){ ## obsolete!!!
        gsea.kegg <- batch.gsea(sort(deg.l[[i]][, 'or']), genesets = hs.gi$kegg)
        gsea.kegg <- data.frame(gsea.kegg, hs.gi$kegg.info[rownames(gsea.kegg), ])
        ## b <- gsea.kegg[, "pvalue"] < 0.05/nrow(gsea.kegg) & gsea.kegg[, "ES"] > 0
        gsea.kegg.l[[i]] <- gsea.kegg
        
        gsea.go <- batch.gsea(sort(deg.l[[i]][, 'or']), genesets = hs.gi$go)
        gsea.go <- data.frame(gsea.go, hs.gi$go.info[rownames(gsea.go), ])
        ##b <- gsea.go[, "pvalue"] < 0.05/nrow(gsea.go) & gsea.go[, "ES"] > 0 & gsea.go[, "ontology"] == 'BP'
        gsea.go.l[[i]] <- gsea.go
      }

      
      if(print.gsea){
        print.n <- 30
        print(hs.gi$kegg.info[rownames(gsea.kegg)[b], , drop=F][1:min(sum(b), print.n), ])
        print(hs.gi$go.info[rownames(gsea.go)[b], 1, drop=F][1:min(sum(b), print.n), ])
      }
    }
    save(deg.l, gsea.go.l, gsea.kegg.l, file=file.path(output.dir, 'middle.res.RData'))
  }
  
  ## load('middle.res.RData')
  
  marker.gene.l <- NULL
  ## my.hclust(raw.d[1:1000, ], sample.class = split(names(hc.10), hc.10))
  for(i in 1:length(deg.l)){
    ts <- deg.l[[i]]
    b.up <- d.grps[rownames(ts), as.character(i)] == apply(d.grps[rownames(ts), ], 1, max)
    b.down <- d.grps[rownames(ts), as.character(i)] == apply(d.grps[rownames(ts), ], 1, min)
    
    if(is.element('or', colnames(ts))){
      ts.up <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (ts[, 'or'] > 2 | is.infinite(ts[,'or'])) & b.up, ]
      ts.down <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (ts[, 'or'] < 0.5 & !is.infinite(ts[,'or'])) & b.down, ]
    }else{
      if(max(ts[, 8], na.rm = T) > 30){
        fc <- log201(ts[, 6]) - log201(ts[, 7])
      }else{
        fc <- ts[,8]
      }
      ts.up <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (fc > 1 & !is.na(fc)) & b.up, ]
      ts.down <- ts[ts[, 'pvalue'] < 0.05/length(ts) & (fc < -1 & !is.na(fc)) & b.down, ]
    }
    
    co.up.gs <- intersect(rownames(gene.info$gene.info), row.names(ts.up))
    co.down.gs <- intersect(rownames(gene.info$gene.info), row.names(ts.down))
    co.gs <- union(co.up.gs, co.down.gs)
    marker.gene.l[[i]] <- data.frame(gene.info$gene.info[co.gs, ,drop=F], ts[co.gs, ,drop=F])
  }
  
  for(i in 1:length(marker.gene.l)){
    co.gs <- rownames(marker.gene.l[[i]])
    if(length(co.gs) == 0){
      next
    }
    td <- norm.d[co.gs, ,drop=F]
    b <- ccv[colnames(td)] == i
    colnames(td)[b] <- '1'
    colnames(td)[!b] <- '.'
    td[is.na(td)] <- 0

    if(nrow(td) > 200){
      td <- td[1:200, ]
    }
    

    width=1200
    height = max(1000 * nrow(td)/2000, 500)
    tiff(file=file.path(output.dir, paste(sampleName, '_heatmap-', i, '.tif', sep = '')), 
           width = width, height = height)
    if(nrow(td) > 1){
      my.heatmap(td, row.clust = F, column.clust = F, 
                 col.center.zero = F, column.label = 'as.is', row.label = NULL,
                 col = color.gradient(low='white', middle = 'blue', high='red', n=100), 
                 columnsep = T, sep.color = 'black',
                 column.class = cell.class.l, X11 = F, main=i)
    }else{
      my.vioplot(split(td[1, ], colnames(td)))  
    }
    dev.off()
  }
  
  all.markers <- sort(unique(unlist(sapply(marker.gene.l, rownames))))
  tiff(file=file.path(output.dir, paste(sampleName, '_heatmap-allMarker.tif', sep='')), 
       height=1000, width = 1200)
  td <- norm.d[all.markers, ,drop=F]
  td[is.na(td)] <- 0
  if(nrow(td) > 0){
    if(nrow(td) > 1){
      row.clust <- TRUE
    }else{
      row.clust <- FALSE
    }
    my.heatmap(td, row.clust = FALSE, column.clust = F, 
               col.center.zero = F, column.label = NULL, 
               col = color.gradient(low='white', middle = 'blue', high='red', n=100), 
               column.class = split(names(ccv), ccv), X11 = F)
  }
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

check.geneset.batch <- function(d, norm.d, gsoi.l, min.sample.n=3, cell.class.vector, output.prefix='this', ...){
  res.l <- list()

  for(gs.name in names(gsoi.l)){
    gs <- gsoi.l[[gs.name]]
    ## x1 <- d[intersect(rownames(d), gs), , drop=F]
    ## b <- apply(x1 > 0, 1, sum) > min.sample.n
    ## raw.d.t <- x1[b, ,drop=F]
    ## norm.d.t <- norm.d[rownames(raw.d.t), , drop=F]
    ## res.l[[gs.name]] <- check.geneset(d = raw.d.t, norm.d = norm.d.t, gs = gs, gs.name =gs.name, 
    res.l[[gs.name]] <- check.geneset(d = d, norm.d = norm.d, gs = gs, gs.name =gs.name, 
                                      cell.class.vector = cell.class.vector, 
                                      output.prefix = paste(output.prefix, gs.name, sep='_'), ...)
  }
  invisible(res.l)
}

## data dir is 10x's data directory containing: barcodes.tsv, genes.tsv, matrix.mtx
read.10x.mtx2Seurat <- function(data.dir, min.genes = 1000, ...){
  library(Seurat)
  t <- Read10X(data.dir = data.dir)
  t2 <- CreateSeuratObject(t, min.genes = min.genes, ...)
  invisible(t2)
}

## d is output of read.10x.mtx
sc.preprocess <- function(d, min.expr.genes = 1000, normalize.method='libsize', dnames=NULL){
  if(is.null(dnames)){
    tt <- get.gene.info(rownames(d))    
    dnames <- tt[, c("alias", "GeneID")]
  }
  d.entrez <- dm.rowname.convert(d, id.m = dnames, keep.dup = TRUE)
  rm(d)
  
  d.entrez <- d.entrez[, apply(d.entrez > 0, 2, sum) > min.expr.genes]
  d.entrez <- d.entrez[, apply(d.entrez, 2, max, na.rm=T) > 10]
  
  mito=check.mito.gene.percentage(d.entrez)
  mito.median <- median(mito$percent.mito)
  mito.mad <- mad(mito$percent.mito)
  mito.cut <- (mito.median + 2* mito.mad)

  if(is.null(mito$percent.mito) || mito.cut == 0){
    rd <- d.entrez
  }else{
    rd <- d.entrez[, names(mito$percent.mito)[mito$percent.mito < mito.cut]]      
  }

  if(is.null(normalize.method)){
    norm.d <- list(method='NULL', nfs = 1, d.norm = rd) 
  }else{
    norm.d <- sc.normalize(rd, method = normalize.method)    
  }
  
  mito.para = list(mito=mito, mito.cut=mito.cut)
  invisible(list(raw.d=rd, norm.d=norm.d, mito.para=mito.para, min.expr.genes=min.expr.genes,
                 gene.alias.info = dnames))
}

## a QC test for classification/clustering results.
test.classes.libsize <- function(raw.d, class.v){
  this.d <- cbind(lib.size=apply(raw.d, 2, sum), class=class.v[colnames(raw.d)])
  res <- anova.simple.wrapper(this.d, y.column = 'lib.size', x.column = 'class', log = F, verbose = T)
  my.vioplot(split(this.d[, 1], this.d[,2]))
  invisible(res)
}

import.10x.to.monocle <- function(data.dir, ...){
  dm <- load_cellranger_matrix(data.dir)
  fd <- fData(dm)
  
  # The number 2 is picked arbitrarily in the line below.
  # Where "2" is placed you should place the column number that corresponds to your
  # featureData's gene short names.
  colnames(fd)[2] <- "gene_short_name"
  cds <- newCellDataSet(exprs(dm),
                            phenoData = new("AnnotatedDataFrame", data = pData(dm)),
                            featureData = new("AnnotatedDataFrame", data = fd),
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())
  invisible(cds)
}

## http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
monocle.format.data <- function(d, si, gi, gi.shortName.column=2, count.data=TRUE){
  library(monocle)
  d[is.na(d)] <- 0
  ## d <- d[apply(d, 1, max) > 0, ]
  si <- si[colnames(d), ]
  gi <- gi[rownames(d), ]
  
  colnames(gi)[gi.shortName.column] <- 'gene_short_name'
  
  pd <- new("AnnotatedDataFrame", data = as.data.frame(si))
  fd <- new("AnnotatedDataFrame", data = as.data.frame(gi))
  
  if(count.data){
    
  }else{
    monocle_data <- new("CellDataSet", exprs = d, 
                        phenoData = pd, 
                        featureData = fd,
                        lowerDetectionLimit = 0.1,
                        expressionFamily = tobit(Lower = 0.1))
    rpc_matrix <- relative2abs(monocle_data, method = "num_genes")
    d <- as(as.matrix(rpc_matrix), "sparseMatrix")
  }
  cds <- new("CellDataSet", exprs = d, 
             phenoData = pd, 
             featureData = fd,
             lowerDetectionLimit = 0.5,
             expressionFamily = negbinomial.size())
  print(-1)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  invisible(cds)
}

monocle.cluster <- function(cds, cluster.method=c("densityPeak", "louvain", "DDRTree")[2], ...){
  stop('do not use this clustering method.')
  cds <- reduceDimension(cds, reduction_method = 'ICA', max_components = 30)
  cds@reducedDimA <- cds@reducedDimS
  cds <- clusterCells(cds, method = cluster.method, k=10, ...)
  plot_cell_clusters(cds)
  invisible(cds)
}


## plot_cell_clusters(cds)
monocle.pseudotime.simple <- function(cds, reduceDimensionMethod = c('ICA', 'DDRTree', 'simplePPT')[2], 
                                      plot.b=TRUE){
  cds <- reduceDimension(cds, max_components = 2, method = reduceDimensionMethod)
  cds <- orderCells(cds)
  plot_cell_trajectory(cds, color_by = "State")
  invisible(cds)
}

## cds.order is output of monocle.pseudotime.simple
monocle.pseudotime.heatmap <- function(cds.order, qval.cutoff=5e-2){
  if(cds.order@dim_reduce_type == 'DDRTree'){
    BEAM_res <- BEAM(cds.order, branch_point = 1, cores = 1)    
  }else{
    BEAM_res <- BEAM(cds.order, cores = 1)
  }

  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

  cds.plot <- cds.order[row.names(subset(BEAM_res, qval < qval.cutoff)),]
  max_cluster = nrow(cds.plot)
  for(n in 2:min(10, max_cluster)){
    
    if(cds.order@dim_reduce_type == 'DDRTree'){
      plot_genes_branched_heatmap(cds.plot, branch_point = 1, num_clusters = n, cores = 1, 
                                  use_gene_short_name = T, show_rownames = T)
    }else{
      plot_genes_branched_heatmap(cds.plot, num_clusters = n, cores = 1, 
                                  use_gene_short_name = T, show_rownames = T)      
    }
  }
  invisible(BEAM_res)
}
## end of Monocle
###############################

sc.gini <- function(norm.d){  ## use to select features to classifying cells
  stop('Do not use this.')
  norm.d[is.na(norm.d)] <- 0
  pos.gini.score <- apply(norm.d, 1, gini)  ## for high expressed markers
  neg.gini.score <- apply(max(norm.d)-norm.d, 1, gini)  ## for low expressed markers
  
  invisible(cbind(pos=pos.gini.score, neg=neg.gini.score))
}

## count.data.deg should be output of sc.deg
batch.gsea.count.data <- function(count.data.deg, value.column='pvalue'){
  if(!exists('hs.gi')){
    load('~/data/geneAnnotation/human/hs.gi.RData')    
  }
  
  cd.deg <- count.data.deg
  b <- (cd.deg[, "group1.gc"] > cd.deg[, "group1.size"]/10 |
          cd.deg[, 'group2.gc'] > cd.deg[, "group2.size"]/10)
  gsea.kegg <- batch.gsea(sort(cd.deg[b, value.column]), genesets = hs.gi$kegg)
  gsea.kegg <- data.frame(gsea.kegg, hs.gi$kegg.info[rownames(gsea.kegg), ])
  
  gsea.go <- batch.gsea(sort(cd.deg[, value.column]), genesets = hs.gi$go)
  gsea.go <- data.frame(gsea.go, hs.gi$go.info[rownames(gsea.go), ])
  
  invisible(list(kegg=gsea.kegg, go.bp=gsea.go[gsea.go[, "ontology"]=='BP', ], 
                 go.cc=gsea.go[gsea.go[, "ontology"]=='CC', ],
                 go.mf=gsea.go[gsea.go[, "ontology"]=='MF', ]))
}

###################
## pseduotime
## 1. Mpath  Nature Communication, 2016
## 2. Monocle
## 3. Wishbone Nature Biotechnology, 2016

########################################################
## analysis after integrated clustering (such as Seurat's anchor-based integration method)
## co <- intersect(names(clusters), names(cells))   ## clusters and cells are vector, with names as cell IDs. 
## xtable <- table(clusters[co], cells[co])  
## group.l <- list(group1=c(samples.1), group2=c(samples.2))
integrated.cluster.followup <- function(xtable, group.l, compare.name=NULL, x11.b=TRUE, ...){
  group1.samples <- group.l[[1]]
  group2.samples <- group.l[[2]]
  group1.name <- names(group.l)[1]
  group2.name <- names(group.l)[2]
  
  if(length(c(group1.samples, group2.samples)) > 1){
    my.heatmap(xtable[, c(group1.samples, group2.samples)], col=color.gradient(low='white'), col.center.zero = F, 
               grid = T, row.label = 'as.is', column.label = 'as.is', column.clust = F, X11 = x11.b, 
               main=paste(compare.name, '::', group1.name, ' vs ', group2.name, sep=''), ...)
  }else{
    plot(0,0)
  }
  
  if(length(group1.samples) > 1 & length(group2.samples) > 1){
    l <- apply(xtable, 1, function(x, ss1, ss2){t.test(x[ss1], x[ss2])}, 
               ss1=group1.samples, ss2=group2.samples)
    res <- sapply(l, function(x) c(x$p.value, x$estimate))
    res <- t(res)
    colnames(res) <- c('p.value', names(group.l))
    res <- res[order(res[, 'p.value']), ]
  }else{
    res <- NA
  }
  invisible(res)
}

sc.edgr <- function(xm, group1.ss, group2.ss, group1.name, group2.name){
  library(edgeR)
  cell2group <- c(rep(group1.name, length(group1.ss)), rep(group2.name, length(group2.ss)))
  names(cell2group) <- c(group1.ss, group2.ss)
  
  xd <- aggregate.II(m = xm[, names(cell2group)], v = cell2group, f = 'sum', byrow = F)
  xd[,1] <- xd[,1]/sum(xd[,1]) * 10^6
  xd[,2] <- xd[,2]/sum(xd[,2]) * 10^6
  d <- DGEList(counts=xm[, names(cell2group)], group=factor(cell2group))
  d <- calcNormFactors(d)
  d1 <- estimateCommonDisp(d, verbose=T)
  res <- exactTest(d1)
  fdr <- p.adjust(res$table$PValue, method = 'fdr')
  res$table <- cbind(res$table, FDR=fdr, xd[rownames(res$table), ])
  res$table <- res$table[order(res$table[, "PValue"]), ]
  res$table <- as.matrix(res$table)
  invisible(res)
}

## wrapper of integrated.cluster.followup()
sc.diff.cell.cluster <- function(x, group1.ss, group2.ss, group1.name, group2.name, 
                                 compare.name=NULL, x11.b=FALSE){
  group1.ss <- intersect(group1.ss, colnames(x))
  group2.ss <- intersect(group2.ss, colnames(x))
  
  two.group.l <- list(group1.ss, group2.ss)
  names(two.group.l) <- c(group1.name, group2.name)
  
  res <- integrated.cluster.followup(x, group.l = two.group.l, x11.b=x11.b, compare.name  = compare.name)
}

read.10x.mtx <- function(matrix.dir, min.genes = 500){
  library(Matrix)
  barcode.path <- file.path(matrix.dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix.dir, "features.tsv.gz")
  matrix.path <- file.path(matrix.dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  mat <- as.matrix(mat)
  if(min.genes > 0){
    b <- apply(mat, 2, function(x) sum(x >0)) >= min.genes
    mat <- mat[, b]
  }
  invisible(mat)
}

check.cluster.identify <- function(cl.marker.l, cell.marker.l, tfs, cluster.d, gi, x11.b=TRUE, ...){
  cl.cl <- my.hclust(log201(cluster.d), main="unsupervised hierarchical clustering of clusters", noX = T)
  cluster.d <- cluster.d[, colnames(cluster.d)[cl.cl$order]]
  
  barplot(sort(sapply(cl.marker.l, length)), main='number of markers in each cluster', 
               ylab='number of marker', xlab='cluster')
  
  cl.gsea.l <- list()
  all.genes <- rownames(cluster.d)
  marker.count.m <- matrix(0, nrow=length(cl.marker.l), ncol=length(cell.marker.l), 
                           dimnames = list(names(cl.marker.l), names(cell.marker.l)))
  for(i in names(cl.marker.l)){
    cl.gsea.l[[i]] <- gsea.fisher(cl.marker.l[[i]], list2matrix(cell.marker.l), all.genes, min.pathway.size = 2, 
                                  alternative = 'greater')
    t <- sapply(cell.marker.l, function(x, y) length(intersect(x,y)), y=cl.marker.l[[i]])
    t <- t/sapply(cell.marker.l, length)
    marker.count.m[i, names(t)] <- t
  }
  
  t <- sapply(cl.gsea.l, function(x) print(x[1:4, ]))
  ## Too many markers appearing in most of clusters, so check the reason:
  x <- sort(table(unlist(cl.marker.l)))
  barplot((table(x)), las=2, xlab='number of appearance in clusters', ylab='number of markers')
  
  if(max(cluster.d) > 20){
    cluster.d <- log201(cluster.d)
  }
  cl.markers <- unique(unlist(cl.marker.l))
  co <- intersect(intersect(cl.markers, rownames(cluster.d)), unlist(cell.marker.l))
  co <- c(co, PTPRC='5788')
  
  if(length(co) > 2){
    my.heatmap(t(cluster.d[co, ]), column.label = gi[co, "Symbol"], row.label = 'as.is', grid = T, 
               col.center.zero = F, col=color.gradient(low='white', middle = 'yellow', high = 'red', n = 100), X11 = x11.b, 
               row.clust = F, main='cell markers', ...)
    my.heatmap(t(cluster.d[co, ]), column.label = gi[co, "Symbol"], row.label = 'as.is', grid = T, 
               col.center.zero = F, col=color.gradient(low='white', 'yellow', high = 'red', n = 100), X11 = x11.b, 
               row.clust = T, main='cell markers', ...)
  }
  
  co <- intersect(intersect(cl.markers, tfs), rownames(cluster.d))
  if(length(co) > 2){
    my.heatmap(t(cluster.d[co, ]), column.label = gi[co, "Symbol"], row.label = 'as.is', grid = T, 
               col.center.zero = F, col=color.gradient(low='white', middle = 'blue', high = 'red', n = 100),
               X11 = x11.b, 
               row.clust = F, main='transcription factors', ...)
    my.heatmap(t(cluster.d[co, ]), column.label = gi[co, "Symbol"], row.label = 'as.is', grid = T, 
               col.center.zero = F, col=color.gradient(low='white', middle = 'blue', high = 'red', n = 100),
               X11 = x11.b, 
               row.clust = T, main='transcription factors', ...)
    
    tfs.marker.m <- matrix(0, nrow=length(cl.marker.l), ncol=length(co), dimnames = list(names(cl.marker.l), co))
    for(tf in co){
      t <- sapply(cl.marker.l, function(x, y) is.element(y, x), y=tf)
      tfs.marker.m[names(t), tf] <- t    
    }
      
    my.heatmap(tfs.marker.m, col=color.gradient('white', middle = 'blue', 'red'), grid = T, 
               row.label = 'as.is', column.label = gi[colnames(tfs.marker.m), "Symbol"], col.center.zero = F, 
               X11 = x11.b, row.clust = F, main='transcription factor appearance')
  }
  
  my.heatmap(marker.count.m > 0, col=color.gradient('white', middle = 'blue', 'red'), grid = T, 
             row.label = 'as.is', column.label = 'as.is', col.center.zero = F, X11 = x11.b, row.clust = F, 
             main='cell marker appearance')
  my.heatmap(marker.count.m > 0, col=color.gradient('white', middle = 'blue', 'red'), grid = T, 
             row.label = 'as.is', column.label = 'as.is', col.center.zero = F, X11 = x11.b, row.clust = T, 
             main='cell marker appearance')
  my.heatmap(marker.count.m, col=color.gradient('white', middle = 'blue', 'red'), grid = T, 
             row.label = 'as.is', column.label = 'as.is', col.center.zero = F, X11 = x11.b, row.clust = T, 
             main='percentage of cell marker appeared')
  
  invisible(list(marker.count.m=marker.count.m, cl.gsea.l=cl.gsea.l))
}

gene2pathway <- function(wd, pathway.l, pathway.size.cutoff = 5){
  b <- sapply(pathway.l, function(x, gs) length(intersect(x, gs)), gs = rownames(wd)) > pathway.size.cutoff 
  pathway.l <- pathway.l[b]
  pm <- sapply(pathway.l, function(x, d){
    co <- intersect(x, rownames(d))
    t2 <- NULL
    if(length(co) > 1){
      t1 <- d[co, ] > 0
      t2 <- apply(t1, 2, function(x, n) sum(x)/n, n=length(co))
    }
    return(t2)
  }, d=wd)
  pm <- t(as.matrix(pm))
  invisible(pm)
}

cluster.pathway.analysis <- function(cl.marker.l, pathway.l, all.genes, min.pathway.size, pathwayName, 
                                     alternative=c("greater", "two.sided", "less")){
  alternative <- match.arg(alternative)
  cl.pathway.l <- list()
  cl.pathway.m <- NULL
  co.pathway <- names(pathway.l)[sapply(pathway.l, length) > min.pathway.size]
  
  for(cluster in names(cl.marker.l)){
    cl.pathway.l[[cluster]] <- gsea.fisher(cl.marker.l[[cluster]], genesets = list2matrix(pathway.l), 
                                           all = rownames(wd), min.pathway.size = min.pathway.size, 
                                           alternative = alternative)
    cl.pathway.m <- cbind(cl.pathway.m, cl.pathway.l[[cluster]][co.pathway, 1])
  }
  colnames(cl.pathway.m) <- names(cl.pathway.l)
  
  write.csv(cl.pathway.m, file=file.path(output.dir, 
                                         paste('cluster', pathwayName, 'enrichment.csv', sep='.')))
  
  invisible(list(l=cl.pathway.l, m=cl.pathway.m))
}


check.geneset.wrapper <- function(wd, gsoi.l, cell.info, comparison.schema.list, 
                                  output.dir, cellType=NULL, title=NULL, max.gene.n=205, ...){
  cell2class <- cell.info[, "cluster"]
  cell2class <- sort(cell2class)
  
  wd <- wd[, names(cell2class)]
  
  gsoi.l <- lapply(gsoi.l, function(x, n) x[1:min(length(x), n)], n=max.gene.n)
  
  check.geneset.batch(d = wd, norm.d = wd, gsoi.l = gsoi.l, min.sample.n = 0, row.clust = T,  
                      cell.class.vector = cell2class, main=paste(title, cellType, 'cluster', sep='_'),
                      output.prefix = file.path(output.dir, paste(cellType, 'cluster_heatmap', sep='_')), 
                      ...)
  
  if(length(unique(cell.info[, "sid"])) > 1){  ##  Single sample analysis, OR NOT?
    for(cmpr in names(comparison.schema.list)){
      pair.name <- paste(names(comparison.schema.list[[cmpr]]), collapse ='_vs_')
      main <- paste(c(title, cellType, cmpr, pair.name), collapse ='_')
      sample2class <- list2vector(comparison.schema.list[[cmpr]])
      cells <-  rownames(cell.info)[is.element(cell.info[, "sid"], names(sample2class))]
      cell2class <- sample2class[cell.info[cells, "sid"]]
      names(cell2class) <- cells
      cell2class <- sort(cell2class)
      
      this.wd <- wd[, names(cell2class)]
      check.geneset.batch(d = this.wd, norm.d = this.wd, gsoi.l = gsoi.l, min.sample.n = 0, row.clust = T,
                          cell.class.vector = cell2class, main=main,
                          output.prefix = file.path(output.dir, paste(main, 'heatmap', sep='_')), ...)
    }
  }
}


sc.single.gene.analysis <- function(wd, umap.d, gene, cell.info, comparison.schema.list, title=NULL, 
                                    cellType=NULL){
  if(is.element(gene, rownames(wd)) & !is.null(title)){
    gid <- gene
  }else{
    gi <- get.gene.info(gene)
    if(nrow(gi) > 1){
      stop('More than 1 record found for the gene.')
    }
    gid <- gi[,"GeneID"]
    gene.name <- gi[, "Symbol"]
    
    if(is.null(title)){
      title <- gene.name
    }
  }
  
  b <- wd[gid, ] > 0
  plot.tsne.II(umap.d, split(names(b), b), main=title)
  
  l <- split(wd[gid, ], cell.info[colnames(wd), "cluster"]) 
  nms <- as.character(sort(unique(cell.info[, "cluster"])))
  l <- l[nms]
  names(l) <- paste('c', names(l), sep='')
  
  my.vioplot(l, main=paste(title, 'cluster', sep='::'), las=2)
  
  means <- sort(sapply(l, mean))
  set.seed(1)
  cols <- sample(rainbow(12), length(means), replace = T)
  barplot(means, col=cols, main=paste(title, 'cluster', sep='::'), las=2)
  
  if(length(unique(cell.info[, "sid"])) > 1){  ##  Single sample analysis, OR NOT?
    for(cmpr in names(comparison.schema.list)){
      pair.name <- paste(names(comparison.schema.list[[cmpr]]), collapse ='_vs_')
      main <- paste(c(title, cellType, cmpr, pair.name), collapse ='_')
      sample2class <- list2vector(comparison.schema.list[[cmpr]])
      cells <-  rownames(cell.info)[is.element(cell.info[, "sid"], names(sample2class))]
      cell2class <- sample2class[cell.info[cells, "sid"]]
      names(cell2class) <- cells
      cell2class <- sort(cell2class)
      
      l <- split(wd[gid, ], cell2class[colnames(wd)])
      
      my.vioplot(l, main=main, las=2)
      
      means <- sort(sapply(l, mean))
      set.seed(1)
      cols <- sample(rainbow(12), length(means), replace = T)
      barplot(means, col=cols, main=main, las=2)
    }
  }
}

plot.umap <- function(umap.d, x, x.range=c(0,1), col=NULL, plot.color.bar=TRUE, ...){
  library(colourschemes)
  if(is.null(col)){
    pal <- color.gradient(low='white', middle = 'blue', high = 'red', n = 100)
    z <- nearestScheme(data.frame(col=pal,values=seq(0,1, length.out = 100)))
    col <- z(x)
  }
  
  cb.x.pos <- ceiling(min(umap.d[, 1])) + 3
  cb.y.pos <- floor(max(umap.d[,2])) - 1
  
  plot(umap.d[names(x), ], col=col, pch='.', ...)
  
  if(plot.color.bar){
    colorbar.plot(cb.x.pos, cb.y.pos, 1:100, strip.length = 0.25, strip.width = 0.05,
                  col=z(seq(x.range[1], x.range[2], length.out = 100)))
  }
}

