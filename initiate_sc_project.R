dirs <- c('Second_analysis', 'First_analysis_LC', 'CD45', 'XCR_bulk_analysis')
files <- c('sample_info.xls', 'CD45/document_CD45.R', 'XCR_bulk_analysis/document_XCR.R', 
           'Second_analysis/document_Second.R')

config.tmpl.dir <- '~/codes/git/SingleCell/'
config.fs <- rbind(c('TEMPLATE_CD45_batch_scRNA_config.R', 'CD45/CD45_batch_scRNA_config.R'), 
                   c('TEMPLATE_XCR_bulk_analysis_config.R', 'XCR_bulk_analysis/XCR_bulk_analysis_config.R'),
                   c('TEMPLATE_Second_batch_scRNA_config.R', 'Second_analysis/Second_batch_scRNA_config.R'), 
                   c('TEMPLATE_document_overall.R', 'document_overall.R')
                   )

for(dr in dirs){
  if(!dir.exists(dr)){
    dir.create(dr)
  }
}

for(f in files){
  if(!file.exists(f)){
    write('', file = f)
  }
}

for(i in 1:nrow(config.fs)){
  if(!file.exists(config.fs[i, 2])){
    file.copy(file.path(config.tmpl.dir, config.fs[i, 1]), config.fs[i, 2])    
  }
}
  
## Done


  