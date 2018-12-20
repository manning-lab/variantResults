library(data.table)
library(ggplot2)

.getChargeRes <- function(){
  bp.cursmk.files <- paste0(res.dir, c("EUR.DBP.CurSmk.2df1.p.001.txt",
                                       "EUR.SBP.CurSmk.2df1.p.001.txt",
                                       "EUR.MAP.CurSmk.2df1.p.001.txt",
                                       "EUR.PP.CurSmk.2df1.p.001.txt"))
  bp.cursmk.files <<- bp.cursmk.files[order(bp.cursmk.files)]
  
  bp.eversmk.files <- paste0(res.dir, c("EUR.SBP.EverSmk.2df1.p.001.txt",
                                        "EUR.MAP.EverSmk.2df1.p.001.txt",
                                        "EUR.PP.EverSmk.2df1.p.001.txt",
                                        "EUR.DBP.EverSmk.2df1.p.001.txt"))
  bp.eversmk.files <<- bp.eversmk.files[order(bp.eversmk.files)]
  
  lipid.cursmk.files <- paste0(res.dir, c("EA.HDL.CurSmk.All.M1.2df1.p.001.txt",
                                          "EA.TG.CurSmk.All.M1.2df1.p.001.txt",
                                          "EA.LDL.CurSmk.All.M1.2df1.p.001.txt"))
  lipid.cursmk.files <<- lipid.cursmk.files[order(lipid.cursmk.files)]
  
  lipid.eversmk.files <- paste0(res.dir, c("EA.HDL.EverSmk.All.M1.2df1.p.001.txt",
                                           "EA.TG.EverSmk.All.M1.2df1.p.001.txt",
                                           "EA.LDL.EverSmk.All.M1.2df1.p.001.txt"))
  lipid.eversmk.files <<- lipid.eversmk.files[order(lipid.eversmk.files)]
  
  all.files <<- c(bp.cursmk.files,bp.eversmk.files,lipid.cursmk.files,lipid.eversmk.files)
}

.getGtexFiles <- function(x){
  gtex.files <<- c("Adipose_Subcutaneous.v7.signif_variant_gene_pairs.txt.gz",
    "Adipose_Visceral_Omentum.v7.signif_variant_gene_pairs.txt.gz",
    "Colon_Sigmoid.v7.signif_variant_gene_pairs.txt.gz",
    "Colon_Transverse.v7.signif_variant_gene_pairs.txt.gz",
    "Liver.v7.signif_variant_gene_pairs.txt.gz",
    "Pancreas.v7.signif_variant_gene_pairs.txt.gz",
    "Small_Intestine_Terminal_Ileum.v7.signif_variant_gene_pairs.txt.gz",
    "Stomach.v7.signif_variant_gene_pairs.txt.gz",
    "Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")
  gtex.names <<- unlist(lapply(gtex.files, function(x) unlist(strsplit(basename(x),"\\."))[1]))
}

makePlot <- function(data, title, pval){
  if (pval == "intP"){
    x_lab <- "Interaction P-value (-log)"
  } else {
    x_lab <- "Trait P-value (-log)"
  }
  epval = "pval_nominal"
  data[,pval] <- -log10(data[,pval])
  data[,epval] <- -log10(data[,epval])
  ymax <- ceiling(max(data[,epval]))
  xmax <- ceiling(max(data[,pval]))
  names(data)[names(data) == pval] <- "pval"
  pval <- "pval"
  plot <- ggplot() +
    geom_point(data = data, aes_string(pval,epval), shape = 1) +
    scale_y_log10(limits = c(4.5,ymax)) +
    xlim(2.5,xmax) +
    ggtitle(title) +
    xlab(x_lab) +
    ylab("eQTL P-value (-log)") +
    scale_shape_discrete(solid=F) +
    geom_hline(yintercept=8, linetype="dashed", color = "blue") +
    geom_vline(xintercept=8, linetype="dashed", color = "red")
  return(plot)
}

mergeGtex <- function(data, gtex, pval = "intP", mark = "MarkerName", log = F){
  data <- merge(data, gtex[,c("MarkerName","gene_id","hgnc_symbol","pval_nominal")], by.x = mark, by.y = "MarkerName")
  if (log == T){
    data$pval_nominal <- -log10(data$pval_nominal)
    data[,pval] <- -log10(data[,pval])
  }
  return(data)
}

loadGTEx <- function(filename, id.type = "ensembl"){
  library(data.table)
  library(biomaRt)
  df <- fread(filename, data.table = F, stringsAsFactors = F)
  chr.pos <- do.call(rbind, lapply(df$variant_id, function(x) unlist(strsplit(x,"_"))[1:2]))
  df$chr <- chr.pos[,1]
  df$chr <- sub("chr","",df$chr)
  df[df$chr == "X", "chr"] <- 23
  df[df$chr == "Y", "chr"] <- 24
  df[df$chr == "M", "chr"] <- 25
  df$chr <- as.numeric(df$chr)
  df$pos <- as.numeric(chr.pos[,2])
  df$gene_id <- unlist(lapply(df$gene_id, function(x) unlist(strsplit(x,"\\."))[1]))
  names(df)[names(df) == "gene_id"] <- "ensembl_gene_id"
  
  mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  if (id.type == "ensembl"){
    coords <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),filters = c("ensembl_gene_id","chromosome_name"), values=list(unique(df$ensembl_gene_id),c(seq(1,22),"X")), mart=mart)  
    df <- merge(df, coords[,c("ensembl_gene_id","hgnc_symbol")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
    names(df)[names(df) == "ensembl_gene_id"] <- "gene_id"
  } else if (id.type == "entrez") {
    coords <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "chromosome_name"),filters = c("ensembl_gene_id","chromosome_name"), values=list(unique(df$ensembl_gene_id),c(seq(1,22),"X")), mart=mart)  
    df <- merge(df, coords[,c("ensembl_gene_id","entrezgene","hgnc_symbol")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
    names(df)[names(df) == "entrezgene"] <- "gene_id"
  } else {
    coords <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),filters = c("ensembl_gene_id","chromosome_name"), values=list(unique(df$ensembl_gene_id),c(seq(1,22),"X")), mart=mart)  
    df <- merge(df, coords[,c("ensembl_gene_id","hgnc_symbol")], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
    names(df)[names(df) == "hgnc_symbol"] <- "gene_id"
  }
  df <- na.omit(df)
  
  df$MarkerName <- paste(df$chr,df$pos, sep = ":")
  return(df)
}


## setup CHARGE summary statistics filepaths
.getChargeRes()

####
.getGtexFiles()

library(ggpubr)
theme_set(theme_pubr())

## load charge res and make plots interaction
titles <- unlist(lapply(all.files, function(x) {paste(unlist(strsplit(basename(x), "\\."))[1:3], collapse = " ")}))

for (i in seq(1,length(gtex.files))){
  print(paste0(i,"/",length(gtex.files)))
  f <- gtex.files[i]
  gtex <- loadGTEx(f)
  plots.gli <- list()
  plots.intp <- list()
  for (j in seq(1,length(all.files))){
    out.file <- sub(".txt",".csv",paste0( gtex.names[i], ".", basename(all.files[j])))
    data <- fread(all.files[j], data.table = F, stringsAsFactors = F)
    data <- mergeGtex(data, gtex)
    fwrite(data, file = out.file, sep = ",", row.names = F)
    plots.gli[[j]] <- makePlot(data, titles[[j]], "P-value")
    plots.intp[[j]] <- makePlot(data, titles[[j]], "intP")
  }
  
  ## make output file
  
  png(file = paste0( gtex.names[i],".v7.gli.dec14.png"),width = 6, height = 84, units = "in", res=200)
  print(
    ggarrange(plotlist = plots.gli, ncol = 1, nrow = 14)
  )
  dev.off()
  
  png(file = paste0(gtex.names[i],".v7.interaction.dec14.png"),width = 6, height = 84, units = "in", res=200)
  print(
    ggarrange(plotlist = plots.intp, ncol = 1, nrow = 14)
  )
  dev.off()
}

