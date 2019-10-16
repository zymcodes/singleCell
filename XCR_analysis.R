args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 1){
  cat('Rscript batch_scRNA.R config_file.R', '\n')
}
config.file <- args[1]
if(!file.exists(config.file)){
  stop('The config file does not exist.')
}

library(diverse)
source('~/codes/affymetrix/affy_survival_procedure.R')
source(config.file)
load(comparison.schema.file)
load(si.file)

output.dir <- file.path(dirname(micxr.res.file), output.dir)
if(!dir.exists(output.dir)){
  dir.create(output.dir)
}

xcr.all.l <- list()
if(length(grep('(RD|RData)$', micxr.res.file)) > 0){
  ne <- new.env()
  x <- load(micxr.res.file, envir = ne)
  xcr.all.l[['All']] <- ne[[x]]
}else{
  for(s in rownames(si)){
    t <- readxl::read_xlsx(micxr.res.file, trim_ws = T, sheet = s)
    t <- as.data.frame(t)
    xcr.all.l[['All']][[s]] <- t
  }
}

for(s in names(xcr.all.l[['All']])){
  t <- xcr.all.l[['All']][[s]]
  if(nrow(t) > 0){
    Vseg <- sapply(strsplit(t[, 'allVHitsWithScore'], split = '*', fixed = T), function(x) x[1])
    Jseg <- sapply(strsplit(t[, 'allJHitsWithScore'], split = '*', fixed = T), function(x) x[1])
    VJ <- paste(Vseg, Jseg, t[, "aaSeqCDR3"], sep=':')
    xcr.all.l[['All']][[s]] <- cbind(t, V=Vseg, J=Jseg, vj=VJ)
  }else{
    xcr.all.l[['All']][[s]] <- cbind(t, V=NULL, J=NULL, vj=NULL)
  }
}

xcr.b.l <- cdr3.len.l <- list()
for(s in names(xcr.all.l[['All']])){
  type.2letter <- substr(xcr.all.l[['All']][[s]][, 'V'], 1, 2)
  type.3letter <- substr(xcr.all.l[['All']][[s]][, 'V'], 1, 3)
  
  type.2s <- sort(unique(type.2letter))
  type.3s <- sort(unique(type.3letter))

  for(type in type.2s){
    type.b <- type.2letter == type
    xcr.b.l[[type]][[s]] <-  type.b
    xcr.all.l[[type]][[s]] <- xcr.all.l[['All']][[s]][type.b, ,drop=F]
    cdr3.len.l[[type]][s] <- list(nchar(unique(xcr.all.l[['All']][[s]][type.b, 'aaSeqCDR3'])))
  }

  for(type in type.3s){
    type.b <- type.3letter == type
    xcr.b.l[[type]][[s]] <-  type.b
    xcr.all.l[[type]][[s]] <- xcr.all.l[['All']][[s]][type.b, , drop=F]
    cdr3.len.l[[type]][s] <- list(nchar(unique(xcr.all.l[['All']][[s]][type.b, 'aaSeqCDR3'])))
  }
  type.b <- type.3letter == 'IGK' | type.3letter == 'IGL'
  xcr.b.l[['IGLight']][[s]] <-  type.b
  xcr.all.l[['IGLight']][[s]] <- xcr.all.l[['All']][[s]][type.b, , drop=F]
  cdr3.len.l[['IGLight']][[s]] <- nchar(unique(xcr.all.l[['All']][[s]][type.b, 'aaSeqCDR3']))
}

pdf(file=file.path(output.dir, 'CDR3_length_distribution.pdf'), width = 12)
for(type in names(xcr.all.l)){
  cdr3s <- unlist(sapply(xcr.all.l[[type]], function(x) unique(x[, "aaSeqCDR3"])))
  cdr3.aa <- unique(cdr3s)
  barplot(table(nchar(cdr3.aa)), main=type, xlab='length of CDR3', ylab='Frequency')
}
dev.off()

cdr3.freq.l <- cdr3.freq.m.l <- list()
pdf(file=file.path(output.dir, 'CDR3_occurency_distribution.pdf'), width=12)
for(type in names(xcr.all.l)){
  for(s in names(xcr.all.l[[type]])){
    x <- xcr.all.l[[type]][[s]]
    l <- sort(sapply(split(x[, 'cloneCount'], x[, "aaSeqCDR3"]), sum), decreasing = T)
    cdr3.freq.l[[type]][s] <- list(l)
    l <- l/sum(l)
    plot(l, type = 'l', xlab='CDR3 ID', ylab='Proportion', main=paste(type, s, sep='::'), pch='*', col=2)
  }
  cdr3.freq.m.l[[type]] <- t(list2matrix.II(cdr3.freq.l[[type]]))
  l.grp <- apply(cdr3.freq.m.l[[type]], 1, sum)
  l.grp <- sort(l.grp/sum(l.grp), decreasing = T)
  plot(l.grp, type = 'l', xlab='CDR3 ID', ylab='Proportion', 
       main=paste(type, 'All samples', sep='::'))
}
dev.off()

## 
pdf(file.path(output.dir, 'XCR_number_in_each_sample.pdf'), width=12)
t <- sapply(xcr.all.l, function(x) sapply(x, nrow))
cm <- list2matrix.II(t)

cmp <- cm[c('All', 'IG', 'TR', 'IGH', 'TRB'), ]
cmp <- cmp[, order(cmp[1, ])]
barplot(cmp, beside = T, col=2:(nrow(cmp)+1))
legend('topleft', legend = rownames(cmp), fill = 2:(nrow(cmp)+1))

cmp <- cm[c('IGH', 'IGK', 'IGL'), ]
cmp <- cmp[, order(cmp[1, ])]
barplot(cmp, beside = T, col=2:(nrow(cmp)+1))
legend('topleft', legend = rownames(cmp), fill = 2:(nrow(cmp)+1))

cmp <- cm[c('TRA', 'TRB', 'TRD', 'TRG'), ]
cmp <- cmp[, order(cmp[1, ])]
barplot(cmp, beside = T, col=2:(nrow(cmp)+1))
legend('topleft', legend = rownames(cmp), fill = 2:(nrow(cmp)+1))

dev.off()

save(xcr.all.l, xcr.b.l, cdr3.len.l, cdr3.freq.l, cdr3.freq.m.l, config.file, 
     si.file, comparison.schema.file, micxr.res.file, output.dir,
     file=file.path(output.dir, 'xcr.preprocessed.RData'))
## Done for preprocessing
#############################

#############################
## subtype anlysis
for(type in names(xcr.all.l)){
  print(type)
  type.output.dir <- file.path(output.dir, type)
  
  if(!dir.exists(type.output.dir)){
    dir.create(type.output.dir)
  }
  
  x.l <- xcr.all.l[[type]]
  
  cdr3.count.allSample <- cdr3.freq.l[[type]]
  
  cdr3.freq.m <- cdr3.freq.m.l[[type]]

  l <- apply(cdr3.freq.m > 0, 1, sum)
  l <- sort(l[l >= 3])
  if(length(l) > 0){
    pdf(file=file.path(type.output.dir, 'hot.CDR3.pdf'), height = max(length(l)/5, 3))
    par(mar=c(5, 15, 4, 2))
    barplot(l, horiz = T, las=2, xlab='Frequency')
    dev.off()
  }
  
  cdr3.freq.test.l <- list()
  enriched.cdr3.l <- list()
  for(cmpr in names(comparison.schema.list)){
    this.class <- list2vector(comparison.schema.list[[cmpr]])
    co.ss <- intersect(colnames(cdr3.freq.m), names(this.class))
    
    if(length(unique(this.class[co.ss]))!=2){
      next()
    }
    
    this.m <- aggregate.II(cdr3.freq.m[, co.ss] > 0, this.class[co.ss], f ='sum', byrow = F)
    grp.length <- table(this.class[co.ss])
    all.length <- length(co.ss)
    
    if(0){
      ps <- NULL
      for(i in 1:nrow(this.m)){
        all.count <- sum(this.m[i, ])
        p.t <- my.phyper.test(this.m[i,1], all.count, all.length - all.count, grp.length[colnames(this.m)[1]])
        ps <- c(ps, min(p.t))
      }
    }else{
      k <- grp.length[colnames(this.m)[1]]
      ps <- apply(this.m, 1, function(x, n, k){min(my.phyper.test(x[1], sum(x), n - sum(x), k))}, 
                  n=all.length, k=k)
    }
    
    this.m <- cbind(this.m, pvalue=ps)
    this.m <- this.m[order(this.m[, 'pvalue']), ]
    cdr3.freq.test.l[[cmpr]] <- this.m
    
    b <- this.m[, "pvalue"] < 0.05
    if(sum(b) > 0){
      enriched.cdr3.l[[cmpr]] <- as.data.frame(this.m[b, ,drop=F])
    }
  }
  
  if(length(enriched.cdr3.l) > 0){
    WriteXLS::WriteXLS(x = enriched.cdr3.l, row.names = T, col.names =T, 
                       ExcelFileName = file.path(type.output.dir, 'enriched.CDR3.xls'))
  }
  
  pdf(file=file.path(type.output.dir, 'Venn_of_enriched.CDR3.pdf'))
  l <- sapply(cdr3.freq.test.l, function(x) rownames(x)[x[, 'pvalue'] < 0.05 & (x[, 1]==0 | x[,2]==0)])
  ll <- sapply(l, length)
  if(sum(ll > 0) <=5 & sum(ll > 0) > 0){
    my.venn(l[ll > 0])  
  }
  dev.clear()
  
  { ## diversity test
    diversity.score.m <- sapply(cdr3.count.allSample, function(x) diversity(data.frame('diversity', names(x), x)))
    pm <- NULL
    for(index in rownames(diversity.score.m)){
      ps <- vector()
      ginis <- unlist(diversity.score.m[index, ])
      for(cmpr in names(comparison.schema.list)){
        co1 <- intersect(names(ginis), comparison.schema.list[[cmpr]][[1]])
        co2 <- intersect(names(ginis), comparison.schema.list[[cmpr]][[2]])
        x1 <- ginis[co1]
        x2 <- ginis[co2]
        x1 <- x1[!is.infinite(x1) & !is.na(x1)]
        x2 <- x2[!is.infinite(x2) & !is.na(x2)]
        if(length(x1) > 1 & length(x2) > 1 & (sd(x1) > 1e-15 | sd(x2) > 1e-15)){ ## diverse::diversity ouput wierd numbers
          ps[cmpr] <- t.test(x1, x2)$p.value          
        }else{
          ps[cmpr] <- NA
        }
      }
      pm <- rbind(pm, ps)
    }
    rownames(pm) <- rownames(diversity.score.m)
    WriteXLS::WriteXLS(as.data.frame(pm), row.names = T, col.names = T, 
                       ExcelFileName = file.path(type.output.dir, 'diversity_t.test.xls'))
    
    ####
    arr.ind <- which(pm < 0.05, arr.ind = T)
    pdf(file=file.path(type.output.dir, 'diversity_comparision_dotBoxplot.pdf'))
    if(nrow(arr.ind) > 0){
      for(i in 1:nrow(arr.ind)){
        index <- rownames(pm)[arr.ind[i, 1]]
        cmpr <- colnames(pm)[arr.ind[i, 2]]
        x <- unlist(diversity.score.m[index, ])
        x.l <- split(x, list2vector(comparison.schema.list[[cmpr]])[names(x)])
        dot.boxplot(x.l, ylab=index, main=cmpr, dot.col = 2:5)
      }
    }
    dev.off()
  }
  
  save(file=file.path(type.output.dir, 'wd.RData'), pm, diversity.score.m, enriched.cdr3.l, 
       cdr3.freq.test.l, arr.ind)
}

cat('\n---------------------\n')
cat('W e l l   d o n e !\n')
cat('---------------------\n\n')

## All done!
##################################

