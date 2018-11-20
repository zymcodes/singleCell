## cytometry data analysis
## https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html
##
## 12 Flow Cytometry Terms And Definitions Most Scientists Get Wrong: https://expertcytometry.com/12-flow-cytometry-terms-and-definitions-most-scientists-get-wrong/
##
## ------------------------------------------------------------
## CellCnn (Arvaniti and Claassen 2016), citrus,  FlowSOM, ConsensusClusterPlus, Rphenograph, destiny, SIMLR, 
##   flowCore, CATALYST.
## 1. Due to its computational requirements, citrus can not be run on entire mass cytometry datasets and one 
##    typically must analyze a subset of the data.
## 2. Neither citrus nor CellCnn are able to directly account for complex designs, such as paired experiments or 
##    presence of batches.
## 3. FlowSOM scales easily to millions of cells and thus no subsetting of the data is required.
## 4. we show how to conduct the differential analysis of cell population abundances using the generalized linear 
##    mixed models (GLMM) and of marker intensities using linear models (LM) and linear mixed models (LMM). Model 
##    fitting is performed with lme4 and stats packages, and hypothesis testing with the multcomp package.
## 5. CATALYST package provides an implementation of the de-barcoding algorithm described by Zunder et al. (Zunder 
##    et al. 2015) and the bead-based normalization from Finck et al. (Finck et al. 2013).
if(0){
  bioc.install('flowCore')
}

read.fcs <- function(fcs.file, ...){
  t <- read.flowSet(fcs.file, transformation = FALSE, truncate_max_range = FALSE, ...)
  expr <- fsApply(t, exprs)
  invisible(expr)
}
