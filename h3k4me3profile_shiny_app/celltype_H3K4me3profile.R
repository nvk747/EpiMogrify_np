
## load libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(grid)
require(gridExtra)
library(data.table)
options(stringsAsFactors = FALSE)


## Merge samples in each celltype. The union of broad peaks across samples are generated.
celltype_unionprofile <- function(peakdir, resdir, celltypes, qValue_cutoff ){
  extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
  celltypes2 <- as.list(read.csv(celltypes,header = FALSE))
  for (ct in celltypes2){
    p = paste0(peakdir,"/", ct)
    files = list.files(p)
    print(files)
    i=0
    if (length(files) != 0){
      x = GRanges()
      cat ("Processing cell type ", ct, "\n" )
      for (f in files){
        cat ("Sample: ", f[1])
        pks = import(paste0(peakdir,"/", ct, "/",f[1]), format = "bed", extraCols = extraCols_broadPeak)
        pks_cutoff = pks[elementMetadata(pks)$qValue  >= qValue_cutoff]
        cat("\t Number of peaks: ", length(pks_cutoff), "\n" )
        if (length(pks_cutoff) >= 10000){
          x = c(x, pks_cutoff) 
        }
      }
      if (length(x) !=0){
        x = reduce(x)
        export.bed(x, paste0(resdir, "/",  ct, "-union.bed"))
      }
    }
  }
}
