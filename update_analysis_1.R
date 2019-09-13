paras <- commandArgs(trailingOnly = TRUE)

if(length(paras) != 2){
  cat("Rscript this.R  processed.res.RData  hclust.res.RData\n")
  quit()
}

processed.data.file <- paras[1]  ## file handle
hclust.data.file <- paras[2]
cat('Rscript this.R ', processed.data.file, hclust.data.file, "\n")

print('Follow and update results of basic_analysis.R.')

source('~/codes/git/SingleCell/SingleCell.R')
source('~/codes/affymetrix/affy_survival_procedure.R')
load('~/data/geneAnnotation/human/ncbi.gene.info.RData')

load(processed.data.file)
load(hclust.data.file)

sampleName <- strsplit(processed.data.file, '_')[[1]][1]
gs.res.l <- list()
load('~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
for(gs.name in names(gsoi.l)){
  gs <- gsoi.l[[gs.name]]
  gs <- gs[is.element(gs, rownames(wd$norm.d$d.norm))]
  if(length(gs) > 200){
    t1 <- apply(wd$norm.d$d.norm[gs, ], 1, max, na.rm=T)
    gs <- gs[is.element(gs, names(sort(t1, decreasing = T))[1:200])]
  }
  ## t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=gs,
  ##                   gs.name = gs.name, cell.class.vector = hcl.10,
  ##                   output.prefix = paste(sampleName, 'gsoi', sep='-'),
  ##                   row.clust = T, display.scale = 'none')
  
  t <- check.geneset(d=wd$raw.d, norm.d = wd$norm.d$d.norm, gs=gs,
                     gs.name = gs.name, cell.class.vector = hcl.10,
                     output.prefix = paste(sampleName, 'gsoi', sep='-'),
                     row.clust = T,  column.clust = T, display.scale = 'none')
  
  gs.res.l[[gs.name]] <- t
}

save(gs.res.l, file=paste(sampleName, '_geneset.res.RData', sep=''))
