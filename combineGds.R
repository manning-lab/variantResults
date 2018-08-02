
library(SeqArray)

# Parse inputs
args <- commandArgs(trailingOnly=T)
gds.files <- unlist(strsplit(args[1], ","))
out.pref <- args[2]

##### test inputs #####
# gds.files <- c("/Users/tmajaria/Documents/projects/mirna.sites.chr22.22.gds")
# out.pref <- "mirna.sites.chr22.all"
#######################


seqMerge(gds.files, paste0(out.pref, ".gds"))
