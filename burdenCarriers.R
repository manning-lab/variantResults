## burdenCarriers.R ##
# R --vanilla --args ${gds_file} ${assoc_file} ${group_file} ${out_pref} ${default="0.001" pval_thresh} < /variantResults/burdenCarriers.R

library(data.table)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(tidyr)

# Parse inputs
args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
assoc.file <- args[2]
group.file <- args[3]
out.pref <- args[4]
pval.thresh <- as.numeric(args[5])

# load association results
assoc.data <- fread(assoc.file, data.table = F, stringsAsFactors = F)

# subset to groups that pass pvalue threshold
assoc.top <- assoc.data[assoc.data[,"pval_0"] < pval.thresh, "MarkerName"]

# load groups (all varaints)
group.data <- fread(group.file, data.table = F, stringsAsFactors = F)

# make sure we have the right name of the column
names(group.data)[names(group.data) == "group.id"] <- "group_id"

# subset to top groups
group.data <- group.data[group.data$group_id %in% assoc.top,]

# open the gds file
gds.data <- seqOpen(gds.file)

# get the current chromosome
chr <- seqGetData(gds.data, "chromosome")[1]

# subset the groups to just this chromosomes groups
group.data <- group.data[group.data$chromosome == chr, ]

# only go on if we have some groups
if (nrow(group.data) > 0) {

  for (gind in seq(1,length(assoc.top))){
    # get the current group
    g <- assoc.top[gind]

    # show me what the current group is
    print(g)

    # load groups
    cur.group <- group.data[group.data$group_id == g,]
    
    # get chromosomes of groups
    min.pos <- min(cur.group$position) - 1
    max.pos <- max(cur.group$position) + 1
    
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
    
    # subset to groups
    seqSetFilterChrom(gds.data, chr, from.bp = min.pos, to.bp = max.pos)
    
    # get variant dataframe to merge for variant id
    var.df <- .expandAlleles(gds.data)
    cur.group <- merge(cur.group, var.df, by.x = c("chromosome", "position","ref","alt"), by.y = c("chromosome", "position","ref","alt"), all.x = T)
    
    # filter by variant id
    seqSetFilter(gds.data, variant.id = cur.group$variant.id)
    
    # get genotypes
    geno <- data.frame(altDosage(gds.data))
    colnames(geno) <- sub("X","v.", colnames(geno))
    geno[is.na(geno)] <- 0
    
    # some refs are wrong so we need to change these allele counts
    mac <- apply(geno,2,sum)
    geno[mac > 2000,] <- 2 - geno[mac > 2000,]
    
    burden <- apply(geno,1,sum)
    
    # get sample ids for carriers
    carriers <- row.names(geno)[burden > 0]
    carriers <- data.frame(sample_id = carriers, burden = burden[burden > 0])

    g <- gsub("/","_",g)

    fwrite(carriers, file = paste(g,"carriers.csv",sep = "."), sep = ",")
    
    # save the dosage matrix
    fwrite(geno, file = paste(g,"genotypes.csv",sep = "."), sep = ",", row.names = T, col.names = T)

    seqClose(gds.data)
  }
} else {
  fwrite(list(),file = paste(g,"genotypes.csv",sep = "."))
  fwrite(list(), file = paste(g,"carriers.csv",sep = "."))
}

