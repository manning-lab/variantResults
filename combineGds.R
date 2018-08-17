# Parse inputs
args <- commandArgs(trailingOnly=T)
gds.files <- unlist(strsplit(args[1], ","))
out.pref <- args[2]

library(SeqArray)
library(SeqVarTools)

##### test inputs #####
gds.files <- c("/Users/tmajaria/Documents/projects/public_workflows/variantResults/test_inputs/mirna.testing.10.gds","/Users/tmajaria/Documents/projects/public_workflows/variantResults/test_inputs/mirna.testing.22.gds")
out.pref <- "mirna.sites.chr10.chr22"
#######################

# this is old and does not work
# seqMerge(gds.files, paste0(out.pref, ".gds"))
#

# need to get all of the data into R data format first
# loop through the input gds files to organize the data
# well only want sample ids, variant ids, positions, chromosomes, alleles, and genotypes

# setup
sample.ids <- c()
variant.ids <- c()
variant.pos <- c()
variant.chr <- c()
variant.allele <- c()
genotypes <- c()


for (gfi in seq(1,length(gds.files))){
  gf <- gds.files[gfi]
  g <- seqOpen(gf)
  s.ids <- seqGetData(g, "sample.id")
  v.ids <- seqGetData(g, "variant.id")
  v.pos <- seqGetData(g, "position")
  v.chr <- seqGetData(g, "chromosome")
  v.allele <- seqGetData(g, "allele")
  g.geno <- data.frame(altDosage(g), stringsAsFactors = F)
  # g.geno$sample <- row.names(g.geno)
  seqClose(g)
  
  sample.ids <- unique(c(sample.ids, s.ids))
  variant.ids <- c(variant.ids, v.ids)
  variant.pos <- c(variant.pos, v.pos)
  variant.chr <- c(variant.chr, v.chr)
  variant.allele <- c(variant.allele, v.allele)
  if (gfi == 1){
    genotypes <- g.geno
  } else {
    genotypes <- merge(genotypes, g.geno, by = "row.names", all = TRUE)
  }
}

# now put things back in the right order
genotypes <- genotypes[match(sample.ids, genotypes$Row.names), paste0("X",variant.ids)]
row.names(genotypes) <- sample.ids
dm <- list()
dm[["sample"]] <- row.names(genotypes)
dm[["variant"]] <- sub("X","",colnames(genotypes))
genotypes <- as.matrix(as.data.frame(lapply(genotypes, as.numeric)))
dimnames(genotypes) <- dm

.repeat_gds <- function(node, elm, count)
{
  val <- rep(elm, 65536L)
  while (count > 0L)
  {
    size <- ifelse(count <= 65536L, count, 65536L)
    append.gdsn(node, val[seq_len(size)], check=FALSE)
    count <- count - size
  }
  invisible()
}

.append_gds <- function(target.node, gdslist, varname)
{
  .MergeNodeAttr(target.node, gdslist, varname)
  for (i in seq_along(gdslist))
    append.gdsn(target.node, index.gdsn(gdslist[[i]], varname))
  readmode.gdsn(target.node)
  invisible()
}


# make a new gds file
gds <- createfn.gds(paste0(out.pref, ".gds"))
on.exit({ closefn.gds(gds) }, add=TRUE)

# get the right header attributes
put.attr.gdsn(gds$root, "FileFormat", "SEQ_ARRAY")
put.attr.gdsn(gds$root, "FileVersion", "v1.0")

# add description
n <- addfolder.gdsn(gds, "description")
put.attr.gdsn(n, "source.format", "SEQARRAY")

# add data
n <- add.gdsn(gds, "sample.id", sample.ids, compress="LZMA_RA",
              closezip=TRUE)
n <- add.gdsn(gds, "variant.id", variant.ids, compress="LZMA_RA",
              closezip=TRUE)
n <- add.gdsn(gds, "position", variant.pos, storage="int32",
              compress="LZMA_RA", closezip=TRUE)
n <- add.gdsn(gds, "chromosome", variant.chr, storage="string",
              compress="LZMA_RA", closezip=TRUE)
n <- add.gdsn(gds, "allele", variant.allele,
              storage="string", compress="LZMA_RA", closezip=TRUE)
n <- addfolder.gdsn(gds, "genotype")
put.attr.gdsn(n, "VariableName", "GT")
put.attr.gdsn(n, "Description", "Genotype")
n1 <- add.gdsn(n, "data", storage="bit2",
               val=genotypes,
               compress="LZMA_RA")
# readmode.gdsn(vg)
n2 <- add.gdsn(n, "@data", storage="uint8",
               compress="LZMA_RA",
               visible=FALSE)
# .repeat_gds(n1, 1L, length(variant.ids))
n2 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
               compress="LZMA_RA", closezip=TRUE)
put.attr.gdsn(n2, "R.colnames",
              c("sample.index", "variant.index", "length"))
add.gdsn(n, "extra", storage="int16", compress="LZMA_RA", closezip=TRUE)

# sync file
sync.gds(gds)

closefn.gds(gds)

gds.data <- seqOpen(paste0(out.pref, ".gds"))
altDosage(gds.data)

nrow(genotypes) # 42 samples
ncol(genotypes) # 32 variants

length(seqGetData(gds.data, "position")) # 32
length(seqGetData(gds.data, "chromosome")) # 32
length(seqGetData(gds.data, "variant.id")) # 32
length(seqGetData(gds.data, "allele")) # 32


length(seqGetData(gds.data, "sample.id")) # 42



# 
# if (bed_flag == 0L)
# {
#   cat(" (transposed)")
#   permdim.gdsn(vg, c(2L,1L))
# }
# .DigestCode(vg, digest, verbose)
# 
# # add a folder for phase information
# if (verbose) cat("    phase")
# n <- addfolder.gdsn(dstfile, "phase")
# 
# n1 <- add.gdsn(n, "data", storage="bit1", valdim=c(nrow(famD), 0L),
#                compress=compress.annotation)
# .repeat_gds(n1, 0L, as.double(nrow(bimD))*nrow(famD))
# readmode.gdsn(n1)
# .DigestCode(n1, digest, TRUE)
# 
# n1 <- add.gdsn(n, "extra.index", storage="int32", valdim=c(3L,0L),
#                compress=compress.annotation, closezip=TRUE)
# put.attr.gdsn(n1, "R.colnames",
#               c("sample.index", "variant.index", "length"))
# add.gdsn(n, "extra", storage="bit1", compress=compress.annotation,
#          closezip=TRUE)
# 
# 
# # add annotation folder
# n <- addfolder.gdsn(dstfile, "annotation")
# 
# # add annotation/id
# n1 <- add.gdsn(n, "id", bimD$snp.id, storage="string",
#                compress=compress.annotation, closezip=TRUE)
# if (verbose) cat("    annotation/id")
# .DigestCode(n1, digest, verbose)
# 
# # add annotation/qual
# n1 <- add.gdsn(n, "qual", storage="float", compress=compress.annotation)
# .repeat_gds(n1, 100.0, nrow(bimD))
# readmode.gdsn(n1)
# if (verbose) cat("    annotation/qual")
# .DigestCode(n1, digest, verbose)
# 
# # add filter
# n1 <- add.gdsn(n, "filter", storage="int32", compress=compress.annotation)
# .repeat_gds(n1, 1L, nrow(bimD))
# readmode.gdsn(n1)
# put.attr.gdsn(n1, "R.class", "factor")
# put.attr.gdsn(n1, "R.levels", c("PASS"))
# put.attr.gdsn(n1, "Description", c("All filters passed"))
# if (verbose) cat("    annotation/filter")
# .DigestCode(n1, digest, verbose)
# 
# # add the INFO field
# addfolder.gdsn(n, "info")
# # add the FORMAT field
# addfolder.gdsn(n, "format")
# 
# # add sample annotation
# if (verbose) cat("    sample.annotation\n")
# n <- addfolder.gdsn(dstfile, "sample.annotation")
# 
# n1 <- add.gdsn(n, "family", famD$FamilyID, compress=compress.annotation,
#                closezip=TRUE)
# .DigestCode(n1, digest, FALSE)
# 
# n1 <- add.gdsn(n, "father", famD$PatID, compress=compress.annotation,
#                closezip=TRUE)
# .DigestCode(n1, digest, FALSE)
# 
# n1 <- add.gdsn(n, "mother", famD$MatID, compress=compress.annotation,
#                closezip=TRUE)
# .DigestCode(n1, digest, FALSE)
# 
# sex <- rep("", length(sample.id))
# sex[famD$Sex==1L] <- "M"; sex[famD$Sex==2L] <- "F"
# n1 <- add.gdsn(n, "sex", sex, compress=compress.annotation, closezip=TRUE)
# .DigestCode(n1, digest, FALSE)
# 
# n1 <- add.gdsn(n, "phenotype", famD$Pheno, compress=compress.annotation,
#                closezip=TRUE)
# .DigestCode(n1, digest, FALSE)
# 
# # sync file
# sync.gds(dstfile)