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
  load(cell.marker.file)
  
  if(exists('keggMeta.file')){
    ne <- new.env()
    x <- load(keggMeta.file, envir = ne)
    keggMeta.l <- ne[[x]]
    load(kegg.file)
  }else{
    keggMeta.l <- NULL
    keggMeta.file <- NULL
    kegg.l <- NULL
    kegg.file <- NULL
  }

  load(gsoi.file)
  umap.d <- as.matrix(read.table(umap.file, as.is=T, sep='\t', header = T, row.names = 1))
  clusters <- as.matrix(read.table(cell.cluster.file, as.is=T, sep='\t', header = T, row.names = 1))[,1]
  cl.marker.m <- as.matrix(read.table(cluster.marker.file, as.is=T, sep='\t', header = T))
  cl.marker.m[,'GeneID'] <- as.character(cl.marker.m[,'GeneID'])
  cl.marker.l <- split(gsub(' ', '', cl.marker.m[, "GeneID"]), cl.marker.m[, "cluster"])
  cl.marker.l <- get.my.geneset.format.batch(cl.marker.l)

  t1 <- t(sapply(strsplit(names(clusters), '_'), function(x) paste(x[-1], collapse = '_')))
  cell.info <- cbind(clusters, t(t1), si[t1, ])
  colnames(cell.info)[1:2] <- c('cluster', 'sid')
  
  ne <- new.env()
  x <- load(data.file, envir = ne)
  wd <-  ne[[x]] ########*****************************************************************
  wd <- as.matrix(wd)
  rm(ne)
  ## colnames(wd) <- gsub('.', '-', colnames(wd), fixed = T)
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

  cell.gene.ns <- apply(wd > 0, 2, sum)
  save(wd, config.file, si, cell.cluster.file, cell.marker.file, keggMeta.file, cell.info, kegg.file, 
       comparison.schema.list, umap.d, output.dir, cell.marker.l, cl.marker.l,  keggMeta.l, kegg.l, 
       cell.gene.ns, 
       file=file.path(output.dir, 'wd.RData'))
}

############
## Analysis:
## 1. bulk analysis:
if(length(unique(cell.info[, "sid"])) > 1){  ##  Single sample analysis, OR NOT?
  xm <- aggregate.II(wd, v=cell.info[colnames(wd), "sid"], f=sum, byrow = F)
  pdf(file=file.path(output.dir, 'bulk.analysis.pdf'), width=15)
  
  my.hclust(log201(xm), main="unsupervised hierarchical clustering of pseudo-bulk samples", noX = T)
  cc <- cor(log201(xm))
  my.heatmap(cc, col.center.zero = F, col = color.gradient(low='white', middle = 'blue', high = 'red'), 
             grid.color = 'grey', grid = T, row.label = 'as.is', column.label = 'as.is', X11 = F)
  dev.off()
  
  edgr.res.l <- list()
  fdr.cutoff <- 0.05
  degs.l <- deg.keggMeta.l <- list()
  for(cmpr in names(comparison.schema.list)){
    this.compare <- comparison.schema.list[[cmpr]]
    edgr.res.l[[cmpr]] <- sc.edgr(xm=xm, group1.ss=intersect(this.compare[[1]], colnames(xm)),
                                  group2.ss = intersect(this.compare[[2]], colnames(xm)),
                                  group1.name = names(this.compare)[1], group2.name = names(this.compare)[2]) 
    deg.b <- edgr.res.l[[cmpr]]$table[, 'FDR'] < fdr.cutoff
    output.m <- edgr.res.l[[cmpr]]$table[deg.b, ]
    output.m <- data.frame(output.m, 
                           gene.info$gene.info[rownames(output.m), 
                                               c("GeneID", "Symbol", "map_location", "description")])
    degs.l[[cmpr]] <- rownames(output.m)
    write.csv(output.m, file=file.path(output.dir, paste('degs', cmpr, 'csv', sep='.')))
    
    if(!is.null(keggMeta.l)){
      deg.keggMeta.l[[cmpr]] <- gsea.fisher(degs.l[[cmpr]], genesets = list2matrix(keggMeta.l), 
                                              all = rownames(wd), min.pathway.size =  10, 
                                              alternative = 'two.sided')
    }
  }
  degs.l <- get.my.geneset.format.batch(degs.l)
  
  deg.keggMeta.m <- t(list2matrix.II(lapply(deg.keggMeta.l, function(x) x[, "fisher.exact.p"])))
  WriteXLS::WriteXLS(x = as.data.frame(deg.keggMeta.m), row.names = T, col.names =T, 
                     ExcelFileName = file.path(output.dir, 'deg.keggMeta.pvalue.matrix.xls'))
  x.l <- lapply(deg.keggMeta.l, as.data.frame)
  WriteXLS::WriteXLS(x = x.l, row.names = T, col.names =T, 
                     ExcelFileName = file.path(output.dir, 'deg.keggMeta.l.xls'))
  
  
  if(length(degs.l) <= 5){
    pdf(file=file.path(output.dir, 'VennDiagram_of_DEGs.pdf'))
    my.venn(degs.l)
    dev.off()
  }

  deg.co.m <- NULL
  
  co <- Reduce(union, sapply(edgr.res.l, function(x, cutoff){ rownames(x$table)[x$table[, 'FDR'] < cutoff]}, 
                                cutoff=fdr.cutoff))
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
  save(config.file, edgr.res.l, co.degs, degs.l, fdr.cutoff, deg.keggMeta.l, deg.keggMeta.m,
        file=file.path(output.dir, 'bulk.degs.RData'))
}else{
  ##  Single sample analysis!
  edgr.res.l <- degs.l <- co.degs <- NULL
}

## 2. anchored clusters and followup analysis
{
  if(!is.null(keggMeta.l)){
    cl.keggMeta.res <- cluster.pathway.analysis(cl.marker.l, pathway.l = keggMeta.l, all.genes = rownames(wd), 
                                         pathwayName = 'keggMeta', min.pathway.size = 10, 
                                         alternative = 'greater')
    cl.kegg.res <- cluster.pathway.analysis(cl.marker.l, pathway.l = kegg.l, all.genes = rownames(wd),
                                      pathwayName = 'kegg', min.pathway.size = 10,
                                      alternative = 'greater')
    wd.kegg <- gene2pathway(wd, kegg.l)
    wd.kegg <- wd.kegg[, rownames(cell.info)[order(cell.info[, "cluster"])]]

    kegg.ps <- intersect(names(kegg.l), rownames(wd.kegg))
    names(kegg.ps) <- kegg.ps
    t1 <- rownames(cl.keggMeta.res$m)[unique(which(cl.keggMeta.res$m < 0.01, arr.ind = T)[,1])]
    ps <- kegg.ps[t1]
    a <- check.geneset(wd.kegg, norm.d = wd.kegg, gs = ps, gs.name = 'kegg metabolism', 
                       plot.each.gene = T, cell.class.vector = cell.info[colnames(wd.kegg), "cluster"], 
                       output.prefix = file.path(output.dir, 'KEGG_metabolism_score'), 
                       main = 'pathway score of kegg metabolism', x11.bulk = F)
    
    t2 <- rownames(cl.kegg.res$m)[unique(which(cl.kegg.res$m < 0.01, arr.ind = T)[,1])]
    ps <- kegg.ps[t2]
    a <- check.geneset(wd.kegg, norm.d = wd.kegg, gs = ps, gs.name = 'kegg', 
                       plot.each.gene = T, cell.class.vector = cell.info[colnames(wd.kegg), "cluster"], 
                       output.prefix = file.path(output.dir, 'KEGG_score'),
                       main = 'pathway score of kegg', x11.bulk = F)
    
    pdf(file=file.path(output.dir, 'KEGG.score_UMAP_sigClusterPathway.pdf'))
    for(pathway in ps[1:min(10, length(ps))]){
      plot.umap(umap.d = umap.d, x=wd.kegg[pathway, ], main=pathway)
    }
    dev.off()
  }
  
  pdf(file=file.path(output.dir, 'cell_gene_number_distribution.pdf'), width = 12)
  ## barplot(sort(table(cell.info[, "group"])), main='groups', ylab='number of cells')
  barplot(sort(table(cell.info[, "cluster"])), main='cluster', ylab='number of cells')
  barplot(sort(table(cell.info[, "sid"])), main='patients', ylab='number of cells', las=2)
 
  mean.genes <- sapply(split(cell.gene.ns, cell.info[names(cell.gene.ns), "cluster"]), mean)
  median.genes <- sapply(split(cell.gene.ns, cell.info[names(cell.gene.ns), "cluster"]), median)
  plot(mean.genes, median.genes[names(mean.genes)],  
       xlab='mean of gene in each cell', ylab='median of genes in each cell', main='gene number')
  text(mean.genes, median.genes[names(mean.genes)], names(mean.genes), adj = 1, cex = 2, col=2)
  
  dev.off()
 
  cluster.d <- aggregate.II(m = wd, v=cell.info[colnames(wd), "cluster"], f = 'sum', byrow = F)
  
  pdf(file=file.path(output.dir, 'cluster.aggregated.analysis.pdf'), width = 15)
  cl.cl <- my.hclust(log201(cluster.d), main="unsupervised hierarchical clustering of clusters", noX = T)
  barplot(sort(table(cell.info[, "cluster"])), las=2)
  dev.off()
  
  #######################################################################
  if(length(unique(cell.info[, "sid"])) > 1){ ##  Single sample analysis, OR NOT?
    x <- table(cell.info[, "cluster"], cell.info[, "sid"])
    x <- t(t(x)/colSums(x))
    
    cl.diff.res.l <- list()
    pdf(file=file.path(output.dir, 'cluster_cell_distribution.pdf'), width=16)
    for(cmpr in names(comparison.schema.list)){
      cl.diff.res.l[[cmpr]] <- integrated.cluster.followup(xtable = x, compare.name = cmpr, 
                                                           group.l = comparison.schema.list[[cmpr]], x11.b = F)
    }
    dev.clear()
    for(cmpr in names(cl.diff.res.l)){
      write.csv(cl.diff.res.l[[cmpr]], 
                file=file.path(output.dir, paste('cluster.diff_', cmpr, '.csv', sep='')))  
    }
    
    degs.20 <- co.degs[, "GeneID"][1:20]
    if(!exists('goi')){
      goi <- degs.20 
    }else{
      goi <- c(goi, degs.20)
    }
    
    cluster.deg.m <- matrix(0, nrow=length(cl.marker.l), ncol=length(degs.l), 
                            dimnames = list(names(cl.marker.l), names(degs.l)))
    for(i in 1:length(cl.marker.l)){
      for(j in 1:length(degs.l)){
        cluster.deg.m[names(cl.marker.l)[i], names(degs.l)[j]] <- length(intersect(cl.marker.l[[i]], degs.l[[j]]))
      }
    }
  }else{
    cluster.deg.m <- goi <- cl.diff.res.l <- NULL
    ##  Single sample analysis!
  }
}


pdf(file=file.path(output.dir, 'cell.identity.pdf'), width=20)
cell.identify.info.l <- check.cluster.identify(cl.marker.l = cl.marker.l, cell.marker.l = cell.marker.l, 
                                               tfs = gsoi.l$tfs, cluster.d = cluster.d, 
                                               gi=gene.info$gene.info, x11.b = FALSE)
dev.off()

save(wd, si, cell.info, config.file, cell.cluster.file, cell.marker.file, keggMeta.file, kegg.file, 
     comparison.schema.list, cell.marker.l, cl.marker.l, keggMeta.l, this.gsoi.l, umap.d, kegg.l,
     cluster.d, cl.diff.res.l, cl.keggMeta.res, cl.kegg.res, cluster.deg.m, cell.identify.info.l, 
     edgr.res.l, degs.l, co.degs, goi, output.dir, wd.kegg,
     file=file.path(output.dir, 'wd.RData'))

## plots
check.geneset.wrapper(wd, this.gsoi.l, cell.info, comparison.schema.list = comparison.schema.list, 
                      output.dir = output.dir, cellType = cellType, x11.bulk=F )
check.geneset.wrapper(wd, degs.l, cell.info, comparison.schema.list = comparison.schema.list, 
                      output.dir = output.dir, cellType = cellType, x11.bulk=F, max.gene.n = 100)

pdf(file=file.path(output.dir, 'gene.analysis.pdf'))
for(g in goi){
  sc.single.gene.analysis(wd, umap.d, g, cell.info, comparison.schema.list)  
}
dev.off()

cat('\n---------------------\n')
cat('W e l l   d o n e !\n')
cat('---------------------\n\n')


## customized analysis
if(0){
  ## data.file <- 'CD19_test/wd.RData'
  ## source(config.file)

  ne <- new.env()
  load(data.file, envir = ne)
  attach(ne)
  custom.gsoi.l <- list()
  check.geneset.wrapper(wd, custom.gsoi.l, cell.info, output.dir = output.dir, cellType = cellType, main = NULL)

  custom.gs <- c()
  pdf(file='custermized.gene.analysis.pdf')
  for(g in custom.gs){
    sc.single.gene.analysis(wd, umap.d, '915', cell.info)    
  }
  detach(ne)
}
