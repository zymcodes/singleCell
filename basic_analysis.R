paras <- commandArgs(trailingOnly = TRUE)

if(length(paras) != 2){
  cat("Rscript this.R  data.dir/data.file  sampleName\n")
  quit()
}
fh <- paras[1]  ## file handle
sampleName <- paras[2]
cat('Rscript this.R ', fh, sampleName, "\n")

source('~/codes/git/SingleCell/SingleCell.R')
source('~/codes/affymetrix/affy_survival_procedure.R')
load('~/data/geneAnnotation/human/ncbi.gene.info.RData')
  
if(file.info(fh)$isdir){  ## if fh is a directory, take the input as 10x genomics output, by default
  rd <- read.10x.mtx(data.dir, min.genes=1000)    
}else{                    ## if fh is a file, take the input as the .xls(matrix) file downloaded from GEO
    if(grepl('RData', fh, ignore.case=T)){
        load(fh)
    }else{    
        rd <- read.table(file = fh, as.is=T, header=T, row.names = 1)
    }
}
print(dim(rd))
wd <- sc.preprocess(rd)
print(dim(wd$norm.d$d.norm))
cell.sizes <- apply(wd$raw.d, 2, sum)
silent.genes <- rownames(wd$raw.d)[apply(wd$raw.d, 1, sum) == 0]
save(wd, cell.sizes, silent.genes, file=paste(sampleName, '_processed.RData', sep=''))

## memory eating !!!
pdf(file=paste(sampleName, '_hclust.pdf', sep=''), width=14)
hcl <- sc.hclust(wd$raw.d, noPlot=FALSE, noX = TRUE, leaflab = "none")
dev.off()
hcl.10 <- cutree(hcl, 10)

hc.res <- sc.hc(d=wd$raw.d, norm.d = wd$norm.d$d.norm, cell.class.vector = hcl.10,
                sampleName = sampleName, gsea.b = F)
save(hc.res, hcl.10, hcl, file=paste(sampleName, '_hclust.RData', sep=''))

## memory eating !!!
tsne.10x.file <- file.path(dirname(dirname(fh)), 'analysis/tsne/2_components/projection.csv')
if(file.exists(tsne.10x.file)){
  tsne.d <- read.table(tsne.10x.file, as.is=T, header=T, row.names = 1, sep=',')
}else{
  ## tsne.d <- tsne::tsne(d)
  d <- wd$norm.d$d.norm[apply(wd$norm.d$d.norm, 1, max, na.rm=T) >0,  ]
  d[is.na(d)] <- 0  
  tsne.res <- Rtsne::Rtsne(t(d), check_duplicates=FALSE, pca=TRUE, dims=2)
  rownames(tsne.res$Y) <- colnames(d)
  tsne.d <- tsne.res$Y
}
pdf(file=paste(sampleName, "_tsne.pdf", sep=''))
plot.tsne(tsne.d, hcl.10[rownames(tsne.d)], fast = FALSE)
dev.off()
## add tsen.d to _hclust.RData
save(hc.res, hcl.10, hcl, tsne.d, tsne.res, file=paste(sampleName, '_hclust.RData', sep=''))


gs.res.l <- list()
load('~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
for(gs.name in names(gsoi.l)){
  gs <- gsoi.l[[gs.name]]
  gs <- gs[is.element(gs, rownames(wd$norm.d$d.norm))]
  if(length(gs) > 200){
    t1 <- apply(wd$norm.d$d.norm[gs, ], 1, max, na.rm=T)
    gs <- gs[is.element(gs, names(sort(t1, decreasing = T))[1:200])]
  }
  t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=gs,
                     gs.name = gs.name, cell.class.vector = hcl.10,
                     output.prefix = paste(sampleName, 'gsoi', sep='-'),
                     row.clust = T, display.scale = 'none')
  
  gs.res.l[[gs.name]] <- t
}

for(gs.name in names(affy.gs)){
  t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=affy.gs[[gs.name]],
                     gs.name = gs.name, cell.class.vector = hcl.10,
                     output.prefix = paste(sampleName, 'affy', sep='-'),
                     row.clust = T, display.scale = 'none')
  
  gs.res.l[[gs.name]] <- t
}

save(gs.res.l, file=paste(sampleName, '_geneset.res.RData', sep=''))
