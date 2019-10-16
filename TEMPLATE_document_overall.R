## Prepare si
if(0){ ## help to make sample_info.xls file.
  load("First_analysis_LC/trans_data/data/exp.Rdata")
  s <- sapply(strsplit(colnames(exp), split='_'), function(x) paste(x[-1], collapse = '_'))
  table(s)
  write.table(cbind(sampleID=sort(unique(s)), patientID='', group=''),sep='\t', row.names = F, 
              file='sample_info.xls')
}


t <- read.table('sample_info.xls', as.is=T, header = T, sep='\t')
si <- as.matrix(t)
rownames(si) <- si[,1]
save(si, file='sample.info.RData')

## Prepare comparison.schema.list
comparison.schema.list <- list(disease=list(cancer=rownames(si)[si[, "group"]=='cancer'], 
                                            normal=rownames(si)[si[, "group"]=='normal']), 
                               location=list(right=rownames(si)[si[, "location"]=='Right'], 
                                             left=rownames(si)[si[, "location"]=='Left']), 
                               cancer.location=list(right=rownames(si)[si[, "location"]=='Right' & si[, "group"]=='cancer'], 
                                                    left=rownames(si)[si[, "location"]=='Left' & si[, "group"]=='cancer']), 
                               normal.location=list(right=rownames(si)[si[, "location"]=='Right' & si[, "group"]=='normal'], 
                                                    left=rownames(si)[si[, "location"]=='Left' & si[, "group"]=='normal'])
                               )
save(file='comparison.schema.list.RData', comparison.schema.list)


## Prepare gsoi.l
load('~/data/geneAnnotation/human/gsoi.l.for.singleCellAnalysis.RData')
this.gsoi.l <- gsoi.l[c("immune", "immune.core", "immuCheck", "immune.cell.marker", 
                        'cell.cycle', 'apoptosis', 'hla', 'dna.repair', 'stem.cell')]
save(file = 'genesets.RData', this.gsoi.l)


## Prepare cell.marker.l
load('~/data/geneAnnotation/human/all.type.of.cell.markers.RData')
save(file='cell.marker.RData', cell.marker.l)

cell.marker.l <- immune.marker.l
save(file='CD45/cell.marker_immune.RData', cell.marker.l)
