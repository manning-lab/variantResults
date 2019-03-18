## burdenCarriers.R ##
# R --vanilla --args ${gds_file} ${variant_file} ${out_pref} ${default="NA" sample_file} ${default="topmedid" id_col} < /variantResults/burdenCarriers.R

# Parse inputs
args <- commandArgs(trailingOnly=T)
gds.files <- unlist(strsplit(args[1],","))
var.file <- args[2]
out.pref <- args[3]
sample.file <- args[4]

library(dplyr)
library(tidyr)
library(data.table)
library(SeqVarTools)

##### test inputs #####
# gds.files <- c("/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr22.pass_and_fail.gtonly.minDP10.gds")
# var.file <- "/Users/tmajaria/Documents/projects/topmed/results/rarevar/freeze5b/noncoding/mirna/by_site/v1/mirna.binding.sites.groups.0.01.small.chr22.csv"
# out.pref <- "mirna.sites.chr22"
# sample.file <- "/Users/tmajaria/Documents/projects/topmed/data/freeze5b_phenotypes/t2d_Model4_cov_age_sex_topmedproject_ALLT2D_clusterAncestryEU.sample.ids.txt"
#######################

# id.col <- args[5]
# pheno.file <- args[4]

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
pos.names <- c("chr","pos","allele","allele.index")
want.names <- c("chromosome","position","alt","minor.allele")
for (i in seq(1,length(pos.names))){
  if (pos.names[i] %in% names(var.data)){
    names(var.data)[names(var.data) == pos.names[i]] <- want.names[i]
  }
}

# stop if we dont have the bare minimum
if (!(all(c("chromosome","position","ref","alt") %in% names(var.data)))){
  stop("Need to have at least chromosome, position, ref, and alt in the variant list file")
}

# check if we have minor allele variable, if not assume alt is minor allele
if (!("minor.allele" %in% names(var.data))){
  var.data$minor.allele <- "alt"
}

# if we have phenotype info, load it
# if (!(pheno.file == "NA")){
#   pheno.data <- fread(pheno.file, data.table = F, stringsAsFactors = F)
# } else {
#   pheno.data <- NA
# }

# if we have sample data, load it
if (!(sample.file == "NA")){
  sample.data <- fread(sample.file, data.table = F, stringsAsFactors = F, header = F)$V1
} else {
  sample.data <- NULL
}

# list to store all carriers
gds.data <- seqOpen(gds.files[1])
all.carriers <- data.frame(sample_id = seqGetData(gds.data, "sample.id"), carrier = 0, stringsAsFactors = F)
seqClose(gds.data)

# list to store all variant ids
all.var.ids <- list()


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
   
    # subset by sample ids if we have them
    if (!is.null(sample.data)){
      seqSetFilter(gds.data, sample.id = sample.data)
    }
      
    # get chromosomes of groups
    min.pos <- min(cur.data$position) - 1
    max.pos <- max(cur.data$position) + 1
      
    # subset to groups
    seqSetFilterChrom(gds.data, chr, from.bp = min.pos, to.bp = max.pos)
    
    # get variant dataframe to merge for variant id
    var.df <- .expandAlleles(gds.data)
    cur.data <- merge(cur.data, var.df, by.x = c("chromosome", "position","ref","alt"), by.y = c("chromosome", "position","ref","alt"), all.x = T)
    
    # filter by variant id
    seqSetFilter(gds.data, variant.id = cur.data$variant.id)
      
    # get genotypes
    geno <- data.frame(altDosage(gds.data))
    colnames(geno) <- sub("X","v.", colnames(geno))
    geno[is.na(geno)] <- 0
      
    # change the dosages for variants if minor allele is ref rather than alt
    var.refs <- cur.data[cur.data$minor.allele == "ref","variant.id"]
    if (length(var.refs) != 0){
      var.refs <- paste0("v.", var.refs)
      for (v in var.refs){
        geno[,v] <- 2 - geno[,v]
      }
    }
    
    burden <- apply(geno,1,sum)
    
    # get sample ids for carriers
    carriers <- row.names(geno)[burden > 0]

    # set carriers to 1
    all.carriers[all.carriers$sample_id %in% carriers, "carrier"] <- 1
    
    # save variant ids so we dont have to look them up later
    all.var.ids[[length(all.var.ids) + 1]] <- cur.data
  }

# close gds
  seqClose(gds.data)
     
}

# combine variant data
all.var.ids <- do.call(rbind,all.var.ids)

# keep only those sample ids that are carriers
all.carriers <- all.carriers[all.carriers$carrier == 1,]

# store gds file paths
new.paths <- c()

# loop back through the genotype files to get a genotype matrix for each chromosome
for (g.f in gds.files){

  # open the gds file
  gds.data <- seqOpen(g.f)

  # set filter to carriers
  seqSetFilter(gds.data, sample.id = all.carriers$sample_id)

  # get the current chromosome
  chr <- seqGetData(gds.data, "chromosome")[1]

  # make sure chr's are in right format
  chr <- as.character(chr)
  cur.data <- all.var.ids[all.var.ids$chromosome == chr,]

  # set filter on gds file
  seqSetFilter(gds.data, variant.id = cur.data$variant.id)#, sample.id = all.carriers$sample_id)
  
  if (length(seqGetData(gds.data, "variant.id")) == 0){
    next
  }

  # define output file name
  out.file <- paste(out.pref, as.character(chr), "gds", sep = ".")
  new.paths <- c(new.paths, out.file)
  
  # write out the gds file
  seqExport(gds.data, out.file)

  # close the gds file
  seqClose(gds.data)

  # clean up
  cleanup.gds(out.file)
  
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