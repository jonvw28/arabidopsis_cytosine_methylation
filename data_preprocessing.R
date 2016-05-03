#
#LIBRARIES
#
source('functions_read_in_first.R')
#
if(!require("plyr")){
        install.packages("plyr")
}
library(plyr)
if(!require("dplyr")){
        install.packages("dplyr")
}
library(dplyr)
#
if(!require("GenomicRanges")){
        source("https://bioconductor.org/biocLite.R")
        biocLite("GenomicRanges")
}
library(GenomicRanges)
if(!require("rtracklayer")){
        source("https://bioconductor.org/biocLite.R")
        biocLite("rtracklayer")
}
library(rtracklayer)
if(!require("stringr")){
        install.packages("stringr")
}
library(stringr)
#
################################################################################
#
#Methylation
#
################################################################################
#
#Load Data
#
load("./Data/At_tiles.RData")
#
# Coerce to data frame and add midpoint
#
At_tiles.data <- Granges.2.df(At_tiles, midpoint = TRUE)
rm(At_tiles)
#
# Define number of columns to keep before inserting new variables
#
breakpoints <- vector("numeric", length = 1)
breakpoints[1] <- 
       as.numeric(length(colnames(dplyr::select(At_tiles.data,seqnames:blast))))
names(breakpoints) <- "data"
#
# Classify according to wild type methylation
#
index <- which.index(At_tiles.data, c("proportionWTCG","proportionWTCHG",
                                      "proportionWTCHH"))

At_tiles.data <- threshold.classes(At_tiles.data, length = breakpoints[1],
                                  index = index, threshold = rep(0.05,times=3),
                                  greater = rep(TRUE, times = 3),
                                  eq = rep(FALSE, times = 3),
                                  head = c("cg_class","chg_class",
                                           "chh_class"))
rm(index)
#
# Keep note of column where new classes end
#
breakpoints <- c(breakpoints,breakpoints[1]+3)
names(breakpoints)[2] <- "classes"
#
#Add Relative Methylation Column for Met1-1 compared to WT
#
index <- which.index(At_tiles.data, c("proportionmet11CG","proportionWTCG"))
At_tiles.data <- col.ratio(At_tiles.data, index = index,
                               split = breakpoints[2], title = "relative_meth")
rm(index)
breakpoints[2] <- breakpoints[2] + 1
#
#class TEL type according to % CpG methylation retained in MET 1-1 mutant vs WT
#returns ETEL = TRUE if retention is below threshold for ETEL, FALSE if above
#returns RTEL = TRUE if retention is above threshold for RTEL, FALSE if below   
#
index <- which.index(At_tiles.data, "relative_meth")
At_tiles.data <- threshold.classes(At_tiles.data, length = breakpoints[2],
                                   index = c(index,index), 
                                   threshold = c(0.05,0.8),
                                   greater = c(FALSE,TRUE),
                                   eq = c(FALSE,FALSE),
                                   head = c("etel_class","rtel_class"))
rm(index)
breakpoints[2] <- breakpoints[2] + 2
#
################################################################################
#
#GENE ANNOTATIONS
#
################################################################################
#
#Load Data
#
gene.annot <- rtracklayer::import("./Data/Arabidopsis_thaliana.TAIR10.27.gff3")
#
#Add column for overlaps with genes
#Add argument type = "within" to get only wholly contatined bins
#
At_tiles.data <- df.overlap(At_tiles.data, gene.annot, split = breakpoints[2],
                            title = "gene_overlap", strand = TRUE)
breakpoints[2] <- breakpoints[2] + 1
temp.name <- names(At_tiles.data)[breakpoints[2] + 1]
#
#Add a tag for each class of gene
#
At_tiles.data <- classes.overlap(At_tiles.data, gene.annot, "type", 
                                 breakpoints[2])
breakpoints <- c(breakpoints, which.index(At_tiles.data,temp.name) - 1)
names(breakpoints)[3] <- "gene_class"
temp.name <- names(At_tiles.data)[breakpoints[3] + 1]
#
#Add tag for biotype
#
At_tiles.data <- classes.overlap(At_tiles.data, gene.annot, "biotype", 
                                 breakpoints[3], append = TRUE, addName = "bio")
breakpoints[3] <- which.index(At_tiles.data,temp.name) -1
#
rm(temp.name,gene.annot)
#
################################################################################
#
#TRANPOSABLE ELEMENT ANNOTATIONS
#
################################################################################
#
#Load data
#
trans.annot <- rtracklayer::import("./Data/Transposable_elements.gff3")
#
#Add column for overlaps wit TE
#Add argument type = "within" to get only wholly contatined bins
#Note need to relabel chromosome names
#
At_tiles.data <- df.overlap(At_tiles.data, trans.annot, split = breakpoints[2],
                            title = "te_overlap", strand = TRUE, relabel = TRUE)
breakpoints[2:3] <- breakpoints[2:3] + 1
temp.name <- names(At_tiles.data)[breakpoints[3] + 1]
#
#
#Add a tag for each class of TE
#
At_tiles.data <- classes.overlap(At_tiles.data, trans.annot, "type", 
                                 breakpoints[3], relabel = TRUE)
breakpoints <- c(breakpoints, which.index(At_tiles.data,temp.name) - 1)
names(breakpoints)[4] <- "te_class"
#
rm(trans.annot,temp.name)
#
################################################################################
#
#EXPRESSION DATA - Where to go next
#
################################################################################
#
#
# Raw data
#
#
files <- c("met13b_a","met13b_b","met13b_c","WT_a","WT_b","WT_c")
pref <- "./Data/"
suff <- "_trimmo_paired_tophat_At_TAIR10_27_sorted_rmdup_picard_genomenorm.bedgraph"
#
for(file in files){
        tmp.exp <- rtracklayer::import(paste(pref,file,suff, sep = ""))
        At_tiles.data <- Granges.score.merge(At_tiles.data,tmp.exp, 
                                             split = breakpoints[2], 
                                             title = paste("exp_",tolower(file),
                                                           sep=""),
                                             cutoff.s = TRUE, cutoff.e = TRUE)
        rm(tmp.exp)
        breakpoints[2:4] <- breakpoints[2:4] + 1
}
rm(files,file,pref,suff)
#
#Now add averages
#
At_tiles.data <- col.gr.mean(At_tiles.data, 
                             index = (breakpoints[2]-5):(breakpoints[2]-3),
                             breakpoints[2] - 3,"avg_exp_met13")
breakpoints[2:4] <- breakpoints[2:4] + 1
At_tiles.data <- col.gr.mean(At_tiles.data, 
                             index = (breakpoints[2]-2):(breakpoints[2]),
                             breakpoints[2],"avg_exp_wt")
breakpoints[2:4] <- breakpoints[2:4] + 1
#
#Also add a relative expression measure
#
At_tiles.data <- col.ratio(At_tiles.data,c(breakpoints[2]-4,breakpoints[2]),
                           breakpoints[2],"rel_exp_met1b_vs_wt")
breakpoints[2:4] <- breakpoints[2:4] + 1
#
################################################################################
#
# Mapability
#
################################################################################
#
files <- c("ath_20mer_0msh","ath_20mer_1msh","ath_20mer_2msh","ath_20mer_3msh",
           "ath_200mer")
pref <- "./Data/"
suff <- ".mappability.bedgraph"
#
for(file in files){
        title <- substring(file, first = 5)
        tmp.map <- rtracklayer::import(paste(pref,file,suff, sep = ""))
        At_tiles.data <- Granges.score.merge(At_tiles.data,tmp.map, 
                                             split = breakpoints[1], 
                                             title = paste("mappab_",title,
                                                           sep=""),
                                             cutoff.s = TRUE, cutoff.e = TRUE)
        rm(tmp.map)
        breakpoints <- breakpoints + 1
}
rm(files,file,pref,suff,title)
#
################################################################################
#
# Annotation
#
################################################################################
#
library(stringr)
#
load("./Data/At_tiles.RData")
gene <- rtracklayer::import("./Data/Arabidopsis_thaliana.TAIR10.27.gff3")
trans <- rtracklayer::import("./Data/Transposable_elements.gff3")

t.df <- as.data.frame(trans)        
regexp <- "[[:digit:]]+"
from <- unique(grep(regexp,t.df[,1],value = TRUE))
to <- stringr::str_extract(from,regexp)
t.df[,1] <- plyr::mapvalues(t.df$seqnames, from = from, to = to)
trans <- GenomicRanges::makeGRangesFromDataFrame(t.df,
                                                 keep.extra.columns = TRUE)
rm(from,to,regexp,t.df)

annotation <- annotation.builder(At_tiles,list(trans,gene,gene),
                                 c(NA,"exon","gene"),
                                 index = c(NA,"type","type"),
                                 names = c("te","exon","intron","intergenic"),
                                 strand = TRUE)
rm(trans,gene,At_tiles)

At_tiles.data <- annot.2.bin(At_tiles.data,annotation,breakpoints[2],
                             strand = TRUE)
breakpoints[2:4] <- breakpoints[2:4] + 1
#
################################################################################
#
saveRDS(At_tiles.data,"data_complete.RData")
