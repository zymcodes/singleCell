## to dos:
## 1. cluster Seurat results, and hierarchically identifye DEGs.

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 1){
  cat('Rscript batch_scRNA.R config_file.R', '\n')
}
config.file <- args[1]
if(!file.exists(config.file)){
  stop('The config file does not exist.')
}

source(config.file)

res.dir <- 'Second_Res'
if(!dir.exists(res.dir)){
  dir.create(res.dir)  
}

output.dir <- file.path(res.dir, output.dir)
if(!dir.exists(output.dir)){
  dir.create(output.dir)
}

this.gene.analysis <- function(wd, umap.d, gene, cell.info, title=NULL){
  gi <- get.gene.info(gene)
  if(nrow(gi) > 1){
    stop('More than 1 record found for the gene.')
  }
  gid <- gi[,"GeneID"]
  gene.name <- gi[, "Symbol"]
  
  if(is.null(title)){
    title <- gene.name
  }
  
  b <- wd[gid, ] > 0
  plot.tsne.II(umap.d, split(names(b), b), main=gene.name)
  
  l <- split(wd[gid, ], cell.info[colnames(wd), "cluster"]) 
  nms <- as.character(sort(unique(cell.info[, "cluster"])))
  l <- l[nms]
  names(l) <- paste('c', names(l), sep='')
  
  my.vioplot(l, main=paste(title, 'cluster', sep='::'), las=2)
  
  means <- sort(sapply(l, mean))
  set.seed(1)
  cols <- sample(rainbow(12), length(means), replace = T)
  barplot(means, col=cols, main=paste(title, 'cluster', sep='::'), las=2)
  
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

this.geneset.analysis <- function(wd, gsoi.l, cell.info, comparison.schema.list, 
                                  output.dir, cellType=NULL, title=NULL, ...){
  cell2class <- cell.info[, "cluster"]
  cell2class <- sort(cell2class)
  
  wd <- wd[, names(cell2class)]
  
  check.geneset.batch(d = wd, norm.d = wd, gsoi.l = gsoi.l, min.sample.n = 0, row.clust = T,  
                      cell.class.vector = cell2class, main=paste(title, cellType, 'cluster', sep='_'),
                      output.prefix = file.path(output.dir, paste(cellType, 'cluster_heatmap', sep='_')), 
                      ...)
  
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


## 0. load codes and metaData
{
  source("~/codes/git/SingleCell/SingleCell.R")
  source("~/codes/affymetrix/affy_survival_procedure.R")
  source('~/codes/pathway/my_gsea.R')
  load('~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
  load('~/data/geneAnnotation/human/ncbi.gene.info.RData')
  load('~/data/pathways/MSigDB/MSigDB_2018.RData')
}

## standard analysis
## input:
## 1. si file; this is global one.
## 2. comparison.schema.list file; this is global one.
## 3. data file
## 4. cell.cluster.file
## 5. umap coordinate file
## 6. genesets of interest file
## 7. output.dir
## 8. cellType name

## 0.1 load data
{
  load(si.file)
  load(comparison.schema.file)
  load(data.file)
  load(gsoi.file)
  umap.d <- as.matrix(read.table(umap.file, as.is=T, sep='\t', header = T, row.names = 1))
  clusters <- as.matrix(read.table(cell.cluster.file, as.is=T, sep='\t', header = T, row.names = 1))[,1]
  
  t1 <- t(sapply(strsplit(names(clusters), '_'), function(x) paste(x[-1], collapse = '_')))
  cell.info <- cbind(clusters, t(t1), si[t1, ])
  colnames(cell.info)[1:2] <- c('cluster', 'sid')
  
  wd <- d.aggr  ########*****************************************************************
  rm(d.aggr)
}

## 0.2 data QC
{
  if(length(setdiff(rownames(cell.info), colnames(wd))) > 0){
    stop("There are some cells in cell.info with no data.")
    if(0){
      cell.info <- cell.info[colnames(wd),]
    }
  }else{
    wd <- wd[, rownames(cell.info)]
  }
  
  if(!setequal(rownames(umap.d), colnames(wd))){
    stop("umap.d and wd have different cells.")
    if(0){
      umap.d <- umap.d[colnames(wd),]
    }
  }
  
  if(length(setdiff(cell.info[, 'sid'], rownames(si))) > 0){
    stop("There are some cells in cell.info do not belong to any sample.")
  }else{
    si <- si[unique(cell.info[, 'sid']), ]
  }
  
  if(length(setdiff(unlist(comparison.schema.list), rownames(si))) > 0){
    warning("There are some samples in comparison.schema.list do not have any annotation.")
    print(setdiff(unlist(comparison.schema.list), rownames(si)))
    for(c1 in names(comparison.schema.list)){
      for(c2 in names(comparison.schema.list[[c1]])){
        comparison.schema.list[[c1]][[c2]] <- intersect(comparison.schema.list[[c1]][[c2]], rownames(si))
      }
    }
  }
}

############
## Analysis:
## 1. bulk analysis:
{
  xm <- aggregate.II(wd, v=cell.info[colnames(wd), "sid"], f=sum, byrow = F)
  pdf(file=file.path(output.dir, 'bulk.analysis.pdf'), width=15)
  
  my.hclust(log201(xm), main="unsupervised hierarchical clustering of pseudo-bulk samples", noX = T)
  cc <- cor(log201(xm))
  my.heatmap(cc, col.center.zero = F, col = color.gradient(low='white', high = 'red'), grid.color = 'grey', 
             grid = T, row.label = 'as.is', column.label = 'as.is', X11 = F)
  dev.off()
  
  edgr.res.l <- list()
  for(cmpr in names(comparison.schema.list)){
    this.compare <- comparison.schema.list[[cmpr]]
    edgr.res.l[[cmpr]] <- sc.edgr(xm=xm, group1.ss=intersect(this.compare[[1]], colnames(xm)),
                                  group2.ss = intersect(this.compare[[2]], colnames(xm)),
                                  group1.name = names(this.compare)[1], group2.name = names(this.compare)[2])    
    deg.b <- edgr.res.l[[cmpr]]$table[, 'PValue'] < 1e-5
    output.m <- edgr.res.l[[cmpr]]$table[deg.b, ]
    output.m <- data.frame(output.m, 
                           gene.info$gene.info[rownames(output.m), 
                                               c("GeneID", "Symbol", "map_location", "description")])
    write.csv(output.m, file=file.path(output.dir, paste('degs', cmpr, 'csv', sep='.')))
  }

  deg.co.m <- NULL
  
  co <- Reduce(union, sapply(edgr.res.l, function(x, cutoff){ rownames(x$table)[x$table[, 'PValue'] < cutoff]}, 
                                cutoff=1e-5))
  t1 <- NULL
  for(cmpr in names(edgr.res.l)){
    if(is.null(t1)){
      t1 <- edgr.res.l[[cmpr]]$table[co, ]
    }else{
      t1 <- cbind(t1, edgr.res.l[[cmpr]]$table[co, ])    
    }
  }
  colnames(t1) <- paste(colnames(t1), rep(names(edgr.res.l), sapply(edgr.res.l, length)),  sep=':')
  
  co.degs <- data.frame(t1, gene.info$gene.info[rownames(t1), c("GeneID", "Symbol", "map_location", "description")])
  write.csv(co.degs, file=file.path(output.dir, 'degs.multi.comparison.COMMON.csv'))
  save(edgr.res.l, co.degs, file=file.path(output.dir, 'bulk.degs.RData'))
}

## 2. anchored clusters and followup analysis
{
  pdf(file=file.path(output.dir, 'cell_number_distribution.pdf'))
  barplot(sort(table(cell.info[, "group"])), main='groups')
  barplot(sort(table(cell.info[, "cluster"])), main='cluster')
  barplot(sort(table(cell.info[, "sid"])), main='patients')
  dev.off()
  
  x <- table(cell.info[, "cluster"], cell.info[, "sid"])
  x <- t(t(x)/colSums(x))
  
  p.res.l <- list()
  pdf(file=file.path(output.dir, 'cluster_cell_distribution.pdf'), width=16)
  for(cmpr in names(comparison.schema.list)){
    p.res.l[[cmpr]] <- integrated.cluster.followup(xtable = x, compare.name = cmpr, 
                                                   group.l = comparison.schema.list[[cmpr]], x11.b = F)
  }
  dev.clear()

  for(cmpr in names(p.res.l)){
    write.csv(p.res.l[[cmpr]], 
              file=file.path(output.dir, paste('cluster.diff_', cmpr, '.csv', sep='')))  
  }
}
  
save(wd, si, cell.cluster.file, cell.info, p.res.l, edgr.res.l, co.degs, comparison.schema.list, umap.d,
     output.dir, file=file.path(output.dir, 'wd.RData'))
  
this.geneset.analysis(wd, this.gsoi.l, cell.info, comparison.schema.list = comparison.schema.list, 
                      output.dir = output.dir, cellType = cellType, x11.bulk=F )

degs.20 <- co.degs[, "GeneID"][1:20]
if(!exists('goi')){
  goi <- degs.20 
}else{
  goi <- c(goi, degs.20)
}

pdf(file=file.path(output.dir, 'gene.analysis.pdf'))
for(g in goi){
  this.gene.analysis(wd, umap.d, g, cell.info)  
}
dev.off()

save(wd, si, cell.cluster.file, cell.info, p.res.l, edgr.res.l, co.degs, comparison.schema.list, umap.d,
     output.dir, goi, file=file.path(output.dir, 'wd.RData'))


## customized analysis
if(0){
  ## data.file <- 'CD19_test/wd.RData'
  ne <- new.env()
  load(data.file, envir = ne)
  attach(ne)
  custom.gsoi.l <- list()
  this.geneset.analysis(wd, custom.gsoi.l, cell.info, output.dir = output.dir, cellType = cellType, main = NULL)

  custom.gs <- c()
  pdf(file='custermized.gene.analysis.pdf')
  for(g in custom.gs){
    this.gene.analysis(wd, umap.d, '915', cell.info)    
  }
  detach(ne)
}
