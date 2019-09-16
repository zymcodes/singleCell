paras <- commandArgs(trailingOnly = TRUE)

if(length(paras) != 4){
  cat("Rscript this.R  data.dir/data.file  sampleName  minGeneNumber geneSet.File\n")
  quit()
}
fh <- paras[1]  ## file handle
sampleName <- paras[2]
min.gene.n <- as.integer(paras[3])
gene.sets.file <- paras[4]
cat('Rscript this.R ', fh, sampleName, min.gene.n, gene.sets.file, "\n")

if(!dir.exists(sampleName)){
  dir.create(sampleName)
}

source('~/codes/git/SingleCell/SingleCell.R')
source('~/codes/affymetrix/affy_survival_procedure.R')
load('~/data/geneAnnotation/human/ncbi.gene.info.RData')

if(file.info(fh)$isdir){  ## if fh is a directory, take the input as 10x genomics output, by default
  rd <- read.10x.mtx(fh, min.genes=min.gene.n)
}else{                    ## if fh is a file, take the input as the .xls(matrix) file downloaded from GEO
  if(grepl('RData', fh, ignore.case=T)){
    load(fh)
  }else{    
    rd <- read.table(file = fh, as.is=T, header=T, row.names = 1)
  }
}
print(dim(rd))

t <- rd[rd > 0][1:10]
xv <- sum(t != ceiling(t))
load('~/data/geneAnnotation/human/ncbi.gene.info.RData')
if(length(intersect(rownames(rd), gene.info$gid2ensg.m[, "Ensembl_ID"])) > min(0.8*nrow(rd), 1000)){
  dnames <- gene.info$gid2ensg.m
}else{
  tt <- get.gene.info(rownames(rd))   
  dnames <- tt[, c("alias", "GeneID")]
}

if(xv > 0){
  print('Not count data.')
  wd <- sc.preprocess(rd, min.expr.genes = min.gene.n, normalize.method = NULL, dnames = dnames)
}else{
  print('Count data.')
  wd <- sc.preprocess(rd, min.expr.genes = min.gene.n, normalize.method = 'libsize', dnames = dnames)
}


print(dim(wd$norm.d$d.norm))
cell.sizes <- apply(wd$raw.d, 2, sum)
silent.genes <- rownames(wd$raw.d)[apply(wd$raw.d, 1, sum) == 0]
save(sampleName, fh, rd, wd, cell.sizes, silent.genes, file=file.path(sampleName, paste(sampleName, '_processed.RData', sep='')))
## the following step take loooooong time!!!
if(xv == 0){
  if(0){  ## too slow. Skip it currently
    wd[['pa.d']] <- get.present.m(wd$raw.d)
    save(sampleName, rd, wd, cell.sizes, silent.genes, file=file.path(sampleName, paste(sampleName, '_processed.RData', sep='')))
  }
}

## memory eating !!!
pdf(file=file.path(sampleName, paste(sampleName, '_hclust.pdf', sep='')), width=14)
## the following step take loooooong time!!!
hcl <- sc.hclust(wd$raw.d, noPlot=FALSE, noX = TRUE, leaflab = "none")
dev.off()
save(hcl, file=file.path(sampleName, paste(sampleName, '_hclust.RData', sep='')))
hcl.10 <- cutree(hcl, 10)
hc.res <- sc.hc(d=wd$raw.d, norm.d = wd$norm.d$d.norm, cell.class.vector = hcl.10,
                sampleName = sampleName, gsea.b = F, output.dir = sampleName)
save(hc.res, hcl.10, hcl, file=file.path(sampleName, paste(sampleName, '_hclust.RData', sep='')))

## memory eating !!!
tsne.10x.file <- file.path(dirname(dirname(fh)), 'analysis/tsne/2_components/projection.csv')
if(file.exists(tsne.10x.file)){
  tsne.d <- read.table(tsne.10x.file, as.is=T, header=T, row.names = 1, sep=',')
}else{
  ## tsne.d <- tsne::tsne(d)
  d <- wd$norm.d$d.norm[apply(wd$norm.d$d.norm, 1, max, na.rm=T) >0,  ]
  d[is.na(d)] <- 0
  library(Rtsne)
  tsne.res <- try(Rtsne(t(d), check_duplicates=FALSE, pca=TRUE, dims=2), FALSE)
  if(length(tsne.res) == 1){
    pca.res <- prcomp(t(d), scale=T, center=T)
    tsne.d <- pca.res$x[, 1:2]
    rownames(tsne.d) <- colnames(d)
    pca.b <- TRUE
  }else{
    rownames(tsne.res$Y) <- colnames(d)
    tsne.d <- tsne.res$Y
    pca.b <- FALSE
  }
}
if(pca.b){
  pdf(file=file.path(sampleName, paste(sampleName, "_PCA.pdf", sep='')))
}else{
  pdf(file=file.path(sampleName, paste(sampleName, "_tsne.pdf", sep='')))
}
for(k in 2:10){
  cls <- cutree(hcl, k)
  plot.tsne(tsne.d, cls[rownames(tsne.d)], fast = FALSE, main=paste(sampleName, k))  
}
dev.off()

## plan to use UMAP replace tsne
if(0){
  d.sr <- CreateSeuratObject(wd$raw.d)
  d.sr@var.genes <- rownames(d.sr@raw.data)
  d.sr <- ScaleData(d.sr)
  d.sr <-RunPCA(d.sr)
  d.sr <- RunUMAP(d.sr)
  umap.d <- d.sr@dr$umap@cell.embeddings
  
  pdf(file=file.path(sampleName, paste(sampleName, "_UMAP.pdf", sep='')))
  for(k in 2:10){
    cls <- cutree(hcl, k)
    plot.tsne(umap.d, cls[rownames(umap.d)], fast=FALSE, main=paste(sampleName, k))
  }
  dev.off()
}


## add tsen.d to _hclust.RData
save(hc.res, hcl.10, hcl, tsne.d, tsne.res, file=file.path(sampleName, paste(sampleName, '_hclust.RData', sep='')))

gs.res.l <- list()
## load('~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
load(gene.sets.file)
for(gs.name in names(gsoi.l)){
  gs <- gsoi.l[[gs.name]]
  gs <- gs[is.element(gs, rownames(wd$norm.d$d.norm))]
  if(length(gs) > 200){
    t1 <- apply(wd$norm.d$d.norm[gs, ], 1, max, na.rm=T)
    gs <- gs[is.element(gs, names(sort(t1, decreasing = T))[1:200])]
  }
  t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=gs,
                     gs.name = gs.name, cell.class.vector = hcl.10,
                     output.prefix = file.path(sampleName, paste(sampleName, 'gsoi', sep='-')),
                     row.clust = T,  column.clust = T, display.scale = 'none', x11.bulk = FALSE)
  
  gs.res.l[[gs.name]] <- t
}

save(gs.res.l, file=file.path(sampleName, paste(sampleName, '_geneset.res.RData', sep='')))

if(0){
  for(gs.name in names(affy.gs)){
    t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=affy.gs[[gs.name]],
                     gs.name = gs.name, cell.class.vector = hcl.10,
                     output.prefix = paste(sampleName, 'affy', sep='-'),
                     row.clust = T, display.scale = 'none')
  
    gs.res.l[[gs.name]] <- t
  }
}
