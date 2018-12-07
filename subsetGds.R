## subsetGds.R ##
# R --vanilla --args ${gds_file} ${variant_file} ${out_pref} < /variantResults/subsetGds.R

# Parse inputs
args <- commandArgs(trailingOnly=T)
gds.files <- unlist(strsplit(args[1],","))
var.file <- args[2]
out.pref <- args[3]

library(dplyr)
library(tidyr)
library(data.table)
library(SeqVarTools)

##### test inputs #####
# gds.files <- c("/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr22.pass_and_fail.gtonly.minDP10.gds")
# var.file <- "/Users/tmajaria/Documents/projects/public_workflows/variantResults/test_inputs/chr22.subsetgds.tsv"
# out.pref <- "/Users/tmajaria/Documents/projects/public_workflows/variantResults/test_outputs/kt2d.chr22.small"
#######################

# functions to generate data frame for variants in group
.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    group_by_("variant.id") %>%
    as.data.frame()
}

# load variants
var.data <- fread(var.file, data.table = F, stringsAsFactors = F)

# make sure we have the right column names
pos.names <- c("chr","pos")
want.names <- c("chromosome","position")
for (i in seq(1,length(pos.names))){
  if (pos.names[i] %in% names(var.data)){
    names(var.data)[names(var.data) == pos.names[i]] <- want.names[i]
  }
}

# stop if we dont have the bare minimum
if (!(all(c("chromosome","position") %in% names(var.data)))){
  stop("Need to have at least chromosome and position in the variant list file")
}

# store gds file paths
new.paths <- c()

# loop through gds files to get all of the sample ids that we'll need
for (g.f in gds.files){

  # open the gds file
  gds.data <- seqOpen(g.f)

  # get the current chromosome
  chr <- seqGetData(gds.data, "chromosome")[1]

  # make sure chr's are in right format
  chr <- as.character(chr)
  var.data$chromosome <- as.character(var.data$chromosome)

  # make sure chrs have same format
  if (startsWith(chr, "chr") & !(startsWith(var.data$chromosome[1], "chr"))){
    var.data$chromosome <- sub("^","chr",var.data$chromosome)
  } else if (!(startsWith(chr, "chr") & startsWith(var.data$chromosome[1], "chr"))){
    var.data$chromosome <- sub("chr","",var.data$chromosome)
  }

  # subset the groups to just this chromosomes groups
  cur.data <- var.data[var.data$chromosome == chr, ]

  
  # only go on if we have some groups
  if (nrow(cur.data) > 0) {
      
    # get chromosomes of groups
    min.pos <- min(cur.data$position) - 1
    max.pos <- max(cur.data$position) + 1
      
    # subset to groups
    seqSetFilterChrom(gds.data, chr, from.bp = min.pos, to.bp = max.pos)
    
    # get variant dataframe to merge for variant id
    var.df <- .expandAlleles(gds.data)
    cur.data <- merge(cur.data, var.df, by.x = c("chromosome", "position"), by.y = c("chromosome", "position"), all.x = T)
    
    # set filter on gds file
    seqSetFilter(gds.data, variant.id = cur.data$variant.id)
    
    if (length(seqGetData(gds.data, "variant.id")) == 0){
      next
    }

    # define output file name
    out.file <- paste(out.pref, as.character(chr), "gds", sep = ".")
    new.paths <- c(new.paths, out.file)
    
    # write out the gds file
    seqExport(gds.data, out.file)
    
    # clean up
    cleanup.gds(out.file)

  }

  # close gds
  seqClose(gds.data)
     
}

if (length(new.paths) > 1){
  # combine gds files
  seqMerge(new.paths, paste0(out.pref, ".gds"), storage.option="LZMA_RA")
} else {
  gds.data <- seqOpen(new.paths)
  seqExport(gds.data, paste0(out.pref, ".gds"))
  seqClose(gds.data)
  cleanup.gds(paste0(out.pref, ".gds"))
}

# delete old files
unlink(new.paths, force = T)
# 