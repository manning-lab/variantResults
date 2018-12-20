############## Test inputs gene eqtl ##############
# genes.toplot <- c("LPL","FADS2","ATRAID,SLC5A6,KRTCAP3,NRBP1")
# gene.file <- "Parker-2017-TableS1.csv"
# data.file <- "EA.TG.CurSmk.All.M1.2df1.txt"
# chr = "chr"
# pos = "pos"
# pval = "P-value"
# gene.id.col = "Ensembl Gene ID"
# gene.id.type = "ensembl"
# out.pref = "EA.TG.CurSmk.Whole.Blood"
# gen = "hg19"
# pad = 10000
# hmm.file = "NA"
# hmm.name = "NA"
# gtex.file <- "Whole_Blood.v7.signif_variant_gene_pairs.txt.gz"
# plotTraitGeneQTL(gene.file, data.file, gtex.file, genes.toplot, chr, pos, pval, gene.id.col, gene.id.type, out.pref, gen, pad, hmm.file, hmm.name)
######################################################

loadChargeRes <- function(filename, pval = "P-value", pt = 1){
  res <- fread(filename, data.table = F, stringsAsFactors = F)
  chr.pos <- do.call(rbind, lapply(res$MarkerName, function(x) unlist(strsplit(x,":"))[1:2]))
  res$chr <- chr.pos[,1]
  res$pos <- as.numeric(chr.pos[,2])
  res[res$chr == "X","chr"] <- 23
  res[res$chr == "Y","chr"] <- 24
  res$chr <- as.numeric(res$chr)
  res <- na.omit(res)
  res <- res[res[,pval] <= pt,]
  return(res)
}

getIntP <- function(data, eff = "IntEffect", std = "IntStdErr"){
  data$intP <- apply(data[,c(eff,std)], 1, function(x) pchisq((x[1]/x[2])^2, df=1, lower.tail=FALSE))
  return(data)
}

loadGenesCoord <- function(gene.file, build = "hg19", id.col = "Entrez.Gene.ID", id.type = "entrez"){
  library(data.table)
  library(biomaRt)
  genes <- fread(gene.file, data.table = F, stringsAsFactors = F)
  if (id.type == "entrez"){
    g.col <- "entrezgene"
  } else if (id.type == "ensembl") {
    g.col <- "ensembl_gene_id"
  } else {
    g.col <- "hgnc_symbol"
  }
  
  if (build == "hg38"){
    mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  } else {
    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }
  
  coords <- getBM(attributes=c("start_position","end_position","chromosome_name", g.col, "hgnc_symbol"),filters = c(g.col,"chromosome_name"), values=list(unique(genes[,id.col]),c(seq(1,22),"X")), mart=mart)  
  if ("Chr" %in% names(genes)){
    genes <- merge(genes, coords[,c("start_position","end_position",g.col,"hgnc_symbol")], by.x = id.col, by.y = g.col)
    names(genes)[names(genes) == "Chr"] <- "chr"
  } else {
    genes <- merge(genes, coords[,c("start_position","end_position","chromosome_name",g.col,"hgnc_symbol")], by.x = id.col, by.y = g.col)
    names(genes)[names(genes) == "chromosome_name"] <- "chr"
  }
  
  return(genes)
}

loadRSCoord <- function(rs.file){
  library(data.table)
  rs <- unique(fread(rs.file, data.table = F, stringsAsFactors = F))
  names(rs)[names(rs) == "SNP.chr"] <- "chr"
  names(rs)[names(rs) == "SNP.Pos"] <- "pos"
  names(rs)[names(rs) == "Inter_Pval"] <- "pval"
  rs[rs$chr == "X","chr"] <- 23
  rs[rs$chr == "Y","chr"] <- 24
  rs[rs$chr == "M","chr"] <- 25
  rs$chr <- as.numeric(rs$chr)
  rs$pos <- as.numeric(rs$pos)
  return(rs)
}

loadRoadmapHMM <- function(filename){
  df <- fread(filename, data.table = F, stringsAsFactors = F)
  df$V1 <- sub("chr","",df$V1)
  df[df$V1 == "X",]$V1 <- 23
  df[df$V1 == "Y",]$V1 <- 24
  df[df$V1 == "M",]$V1 <- 25
  df <- df[,c(1,2,3,4,9)]
  names(df) <- c("chr","start","end","state","color")
  
  # rgb to hex
  s.c <- unique(df[,c("state","color")])
  for (i in seq(1,nrow(s.c))){
    vec <- as.numeric(unlist(strsplit(s.c$color[i],",")))
    s.c[i,"color"] <- rgb(vec[1],vec[2],vec[3], maxColorValue = 255)
  }
  
  df <- merge(df[,c("chr","start","end","state")], s.c, by.x = "state", by.y = "state")
  df$chr <- as.numeric(df$chr)
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  return(df)
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
  return(df)
}

plotTraitGeneOnly <- function(gene.file, data.file, genes.toplot, chr = "chr", pos = "pos", pval = "P-value", gene.id.col = "Entrez.Gene.ID", gene.id.type = "entrez", out.pref, gen = "hg19", pad = 500000, hmm.file = "NA", hmm.name = "ChromHMM"){
  
  # load my required library
  library(data.table)
  library(GenomicRanges)
  library(Gviz)
  library(biomaRt)
  
  # load genes
  genes <- loadGenesCoord(gene.file, build = gen, id.col = gene.id.col, id.type = gene.id.type)
  
  # if hhm, load
  if (hmm.file != "NA"){
    hmm <- loadRoadmapHMM(hmm.file)  
  } else {
    hmm <- c()
  }
  
  # load data
  all.data <- loadChargeRes(data.file, pval = pval)
  all.data <- getIntP(all.data)
  
  # mart
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
  # loop through
  for (gene.names.tosplit in genes.toplot){
    gene.names <- unlist(strsplit(gene.names.tosplit,","))
    
    # get positions for genes of interest
    gene.coords <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),filters = c("hgnc_symbol","chromosome_name"), values=list(gene.names,c(seq(1,22),"X")), mart=mart)  
    
    # get start and end positions
    min.gene.pos <- min(gene.coords$start_position) - 500000
    max.gene.pos <- max(gene.coords$end_position) + 500000
    
    min.pos <- min.gene.pos
    max.pos <- max.gene.pos
    
    # subset to variants in group
    data <- all.data[all.data[,chr] == gene.coords$chr[1] & all.data[,pos] >= min.pos & all.data[,pos] <= max.pos,]
    
    # subset p-values
    data[data[,pval] < 1e-20, pval] <- 1e-20
    data[data[,"intP"] < 1e-20, "intP"] <- 1e-20
    
    # y limit 
    ymin <- floor(min(min(-log10(data[,pval])), min(-log10(data[,"intP"]))))
    ymax <- floor(max(max(-log10(data[,pval])), max(-log10(data[,"intP"]))))
    ylim <- c(ymin, ymax+1)
    
    # make gli track
    gli.track <- DataTrack(
      data = -log10(data[,pval]), 
      start = data[,pos], 
      end = data[,pos], 
      chromosome = data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Trait -log10(p)", 
      background.title="orangered4",
      col = "orangered4",
      fill = "orangered4",
      col.frame="orangered4",
      frame = TRUE
    )
    
    # make interaction track
    int.track <- DataTrack(
      data = -log10(data[,"intP"]), 
      start = data[,pos], 
      end = data[,pos], 
      chromosome = data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Interaction -log10(p)", 
      background.title="cadetblue4",
      col = "cadetblue4",
      fill = "cadetblue4",
      col.frame="cadetblue4",
      frame = TRUE
    )
    
    # gene track
    gtrack <- BiomartGeneRegionTrack(
      start = min.pos, 
      end = max.pos, 
      biomart = mart, 
      strand = "+-", 
      genome = gen, 
      chromosome = data[1,chr], 
      name = "Genomic context", 
      showId = TRUE, 
      geneSymbol = FALSE,
      collapseTranscripts = TRUE,
      background.title = "chartreuse4"
    )
    
    # highlight the current gene
    highlight.track <- HighlightTrack(
      trackList = list(gtrack),
      start = genes[genes$chr == data[1,chr],]$start_position, 
      end = genes[genes$chr == data[1,chr],]$end_position,
      chromosome = data[1,chr]
    )
    
    # hmm track
    if (hmm.file != "NA"){
      # subset hmm
      cur.hmm <- hmm[hmm$chr == cur.data[1,chr] & hmm$end >= min.pos & hmm$start <= max.pos,]
      hmm.track <- AnnotationTrack(
        name = hmm.name, 
        start = cur.hmm$start, 
        end = cur.hmm$end, 
        id = cur.hmm$state,
        feature = cur.hmm$color,
        chromosome = cur.hmm$chr, 
        genome = gen,
        stacking="dense",
        col="transparent",
        groupAnnotation="feature",
        "#FFFFFF"="#FFFFFF",
        "#FF0000"="#FF0000",
        "#CD5C5C"="#CD5C5C",
        "#FFC34D"="#FFC34D",
        "#BDB76B"="#BDB76B",
        "#C0C0C0"="#C0C0C0",
        "#E9967A"="#E9967A",
        "#808080"="#808080",
        "#FFFF00"="#FFFF00",
        "#8A91D0"="#8A91D0",
        "#006400"="#006400",
        "#FF4500"="#FF4500",
        "#008000"="#008000",
        "#C2E105"="#C2E105",
        "#66CDAA"="#66CDAA",
        "#32CD32"="#32CD32"
      )
    }
    
    # scale track
    axTrack <- GenomeAxisTrack(genome = gen , chromosome = data[1,chr])
    
    # idiogram track
    idxTrack <- IdeogramTrack(genome = gen , chromosome = data[1,chr])
    
    # outfile
    out.file <- paste(out.pref, paste(gene.names, collapse = "."), ".png", sep = ".")
    
    # open png
    png(out.file, height = 2*length(gene.names)+5, width = 11, units = "in", res = 400)
    
    # plot it all
    if (hmm.file == "NA"){
      plotTracks(
        c(idxTrack, axTrack, gli.track, int.track, highlight.track),
        showTitle = TRUE,
        sizes=c(1, 2, 5, 5, 5),
        from = min.pos,
        to = max.pos
      )
    } else {
      plotTracks(
        c(idxTrack, axTrack, gli.track, int.track, highlight.track, hmm.track),
        showTitle = TRUE,
        sizes=c(1, 2, 5, 5, 5, 1),
        from = min.pos,
        to = max.pos
      )
    }
    dev.off()
  }
}


plotTraitGeneQTL <- function(gene.file, data.file, gtex.file, genes.toplot, chr = "chr", pos = "pos", pval = "P-value", gene.id.col = "Entrez.Gene.ID", gene.id.type = "entrez", out.pref, gen = "hg19", pad = 500000, hmm.file = "NA", hmm.name = "ChromHMM"){
  
  # load my required library
  library(data.table)
  library(GenomicRanges)
  library(Gviz)
  library(biomaRt)
  
  # load genes
  genes <- loadGenesCoord(gene.file, build = gen, id.col = gene.id.col, id.type = gene.id.type)
  
  # if hhm, load
  if (hmm.file != "NA"){
    hmm <- loadRoadmapHMM(hmm.file)  
  } else {
    hmm <- c()
  }
  
  # load gtex data
  gtex.full <- loadGTEx(gtex.file, gene.id.type)
  
  # load data
  all.data <- loadChargeRes(data.file, pval = pval)
  all.data <- getIntP(all.data)
  
  # mart
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
  # loop through
  for (gene.names.tosplit in genes.toplot){
    gene.names <- unlist(strsplit(gene.names.tosplit,","))
    
    # subset gtex to those that we want
    gtex <- gtex.full[gtex.full$hgnc_symbol %in% gene.names,]
    
    # get positions for genes of interest
    gene.coords <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),filters = c("hgnc_symbol","chromosome_name"), values=list(gene.names,c(seq(1,22),"X")), mart=mart)  
    
    # get start and end positions
    min.gene.pos <- min(gene.coords$start_position) - 500000
    max.gene.pos <- max(gene.coords$end_position) + 500000
    min.gtex.pos <- min(gtex$pos) - 1000
    max.gtex.pos <- max(gtex$pos) + 1000
    min.pos <- min(min.gene.pos, min.gtex.pos)
    max.pos <- max(max.gene.pos, max.gtex.pos)
    
    # subset to variants in group
    data <- all.data[all.data[,chr] == gene.coords$chr[1] & all.data[,pos] >= min.pos & all.data[,pos] <= max.pos,]
    
    # subset p-values
    data[data[,pval] < 1e-20, pval] <- 1e-20
    data[data[,"intP"] < 1e-20, "intP"] <- 1e-20
    gtex[gtex$pval_nominal < 1e-20, "pval_nominal"] <- 1e-20
    
    # y limit 
    ymin <- floor(min(min(-log10(data[,pval])), min(-log10(gtex$pval_nominal)), min(-log10(data[,"intP"]))))
    ymax <- floor(max(max(-log10(data[,pval])), max(-log10(gtex$pval_nominal)), max(-log10(data[,"intP"]))))
    ylim <- c(ymin, ymax+1)
    
    # make gtex tracks
    gtex.tracks <- list()
    cols <- c("black","red","darkgreen","blue")
    i <- 1
    for (g in gene.names){
      g.tex <- gtex[gtex$hgnc_symbol == g,]
      gtex.tracks[[i]] <- DataTrack(
        data = -log10(g.tex$pval_nominal), 
        start = g.tex[,pos], 
        end = g.tex[,pos], 
        chromosome = g.tex[1,chr], 
        genome = gen, 
        ylim = ylim,
        name = paste(g, "eQTL -log10(p)", sep = " "), 
        background.title="black",
        col = cols[i],
        pch = 1,
        col.frame=cols[i],
        frame = TRUE
      )
      i <- i+1
    }
    
    # make gli track
    gli.track <- DataTrack(
      data = -log10(data[,pval]), 
      start = data[,pos], 
      end = data[,pos], 
      chromosome = data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Trait -log10(p)", 
      background.title="orangered4",
      col = "orangered4",
      fill = "orangered4",
      col.frame="orangered4",
      frame = TRUE
    )
    
    # make interaction track
    int.track <- DataTrack(
      data = -log10(data[,"intP"]), 
      start = data[,pos], 
      end = data[,pos], 
      chromosome = data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Interaction -log10(p)", 
      background.title="cadetblue4",
      col = "cadetblue4",
      fill = "cadetblue4",
      col.frame="cadetblue4",
      frame = TRUE
    )
    
    # bars connecting gtex and summary variants
    snp.gtex.track <- HighlightTrack(
      trackList = c(gtex.tracks, list(gli.track, int.track)),
      start = gtex[gtex$pval_nominal < 1e-10,]$pos - 2,
      end = gtex[gtex$pval_nominal < 1e-10,]$pos + 2,
      chromosome = data[1,chr],
      col="lightgrey",
      fill="lightgrey",
      alpha=0.7
    )
    
    # gene track
    gtrack <- BiomartGeneRegionTrack(
      start = min.pos, 
      end = max.pos, 
      biomart = mart, 
      strand = "+-", 
      genome = gen, 
      chromosome = data[1,chr], 
      name = "Genomic context", 
      showId = TRUE, 
      geneSymbol = FALSE,
      collapseTranscripts = TRUE,
      background.title = "chartreuse4"
    )
    
    # highlight the current gene
    highlight.track <- HighlightTrack(
      trackList = list(gtrack),
      start = genes[genes$chr == data[1,chr],]$start_position, 
      end = genes[genes$chr == data[1,chr],]$end_position,
      chromosome = data[1,chr]
    )
    
    # hmm track
    if (hmm.file != "NA"){
      # subset hmm
      cur.hmm <- hmm[hmm$chr == cur.data[1,chr] & hmm$end >= min.pos & hmm$start <= max.pos,]
      hmm.track <- AnnotationTrack(
        name = hmm.name, 
        start = cur.hmm$start, 
        end = cur.hmm$end, 
        id = cur.hmm$state,
        feature = cur.hmm$color,
        chromosome = cur.hmm$chr, 
        genome = gen,
        stacking="dense",
        col="transparent",
        groupAnnotation="feature",
        "#FFFFFF"="#FFFFFF",
        "#FF0000"="#FF0000",
        "#CD5C5C"="#CD5C5C",
        "#FFC34D"="#FFC34D",
        "#BDB76B"="#BDB76B",
        "#C0C0C0"="#C0C0C0",
        "#E9967A"="#E9967A",
        "#808080"="#808080",
        "#FFFF00"="#FFFF00",
        "#8A91D0"="#8A91D0",
        "#006400"="#006400",
        "#FF4500"="#FF4500",
        "#008000"="#008000",
        "#C2E105"="#C2E105",
        "#66CDAA"="#66CDAA",
        "#32CD32"="#32CD32"
      )
    }
    
    # scale track
    axTrack <- GenomeAxisTrack(genome = gen , chromosome = data[1,chr])
    
    # idiogram track
    idxTrack <- IdeogramTrack(genome = gen , chromosome = data[1,chr])
    
    # outfile
    out.file <- paste(out.pref, paste(gene.names, collapse = "."), ".png", sep = ".")
    
    # open png
    png(out.file, height = 2*length(gene.names)+5, width = 11, units = "in", res = 400)
    
    # plot it all
    if (hmm.file == "NA"){
      plotTracks(
        c(idxTrack, axTrack, snp.gtex.track, highlight.track),
        showTitle = TRUE,
        sizes=c(1, 2, rep(5,length(gene.names)), 5, 5, 5),
        from = min.pos,
        to = max.pos
      )
    } else {
      plotTracks(
        c(idxTrack, axTrack, snp.gtex.track, highlight.track, hmm.track),
        showTitle = TRUE,
        sizes=c(1, 2, rep(5,length(gene.names)), 5, 5, 5, 1),
        from = min.pos,
        to = max.pos
      )
    }
    dev.off()
  }
}


plotTraitGeneGtex <- function(gene.file, data.file, gtex.file, pthresh = 1, chr = "chr", pos = "pos", pval = "P-value", gene.id.col = "Entrez.Gene.ID", gene.id.type = "entrez", out.pref = "", gen = "hg19", pad = 500000, hmm.file = "NA", hmm.name = "ChromHMM"){
  
  # load my required library
  library(data.table)
  library(GenomicRanges)
  library(Gviz)
  library(biomaRt)
  
  # load genes
  genes <- loadGenesCoord(gene.file, build = gen, id.col = gene.id.col, id.type = gene.id.type)
  
  # if hhm, load
  if (hmm.file != "NA"){
    hmm <- loadRoadmapHMM(hmm.file)  
  } else {
    hmm <- c()
  }
  
  # load gtex data
  gtex <- loadGTEx(gtex.file, gene.id.type)
  
  # load variant data
  data <- loadChargeRes(data.file, pval = pval, pt = pthresh)
  
  # counter
  i <- 1
  num.genes <- length(unique(genes[,gene.id.col]))
  
  # loop through and plot each desired group
  for (g in unique(genes[,gene.id.col])){
    print(paste(as.character(i), as.character(num.genes), sep = "/"))
    i <- i+1
    
    # subset to just this gene
    g.dat <- genes[genes[,gene.id.col] == g,][1,]
    
    # output file name
    if ("hgnc_symbol" %in% names(genes)){
      g.name <- g.dat$hgnc_symbol
    } else {
      g.name <- g
    }
    
    out.file <- paste(out.pref, g.name, "png", sep = ".")  
    
    # subset to variants in group
    cur.data <- data[data[,chr] == g.dat$chr & data[,pos] >= (g.dat$start_position - pad) & data[,pos] <= (g.dat$end_position + pad), ]
    
    if (nrow(cur.data) == 0) next
    if (nrow(cur.data[cur.data[,pval] < 10e-8,]) == 0) next
    
    # subset gtex to gene
    cur.gtex <- gtex[gtex$gene_id == g,]
    
    if (nrow(cur.gtex) == 0){
      # set out min and max positions to be updated as we go
      min.pos <- g.dat$start_position - pad
      max.pos <- g.dat$end_position + pad
      
      # y limit 
      ylim <- c(floor(min(-log10(cur.data[,pval]), na.rm = T)), ceiling(max(-log10(cur.data[,pval]), na.rm = T)))
    } else {
      # set out min and max positions to be updated as we go
      min.pos <- min(g.dat$start_position - pad, cur.gtex$pos - 2000)
      max.pos <- max(g.dat$end_position + pad, cur.gtex$pos + 2000)
      
      # y limit 
      ylim <- c(floor( min(min(-log10(cur.data[,pval]), na.rm = T),min(-log10(cur.gtex$pval_nominal), na.rm = T) ) ), ceiling( max(max(-log10(cur.data[,pval]), na.rm = T), max(-log10(cur.gtex$pval_nominal), na.rm = T))))
      
      # make snp track
      gtex.track <- DataTrack(
        data = -log10(cur.gtex$pval_nominal), 
        start = cur.gtex[,pos], 
        end = cur.gtex[,pos], 
        chromosome = cur.data[1,chr], 
        genome = gen, 
        ylim = ylim,
        name = paste(g.name, "eQTL -log10(p)", sep = " "), 
        background.title="blue",
        col = "blue",
        pch = 17,
        col.frame="blue",
        frame = TRUE
      )
      
    }
    
    # fix boundaries for snp track
    if (min.pos < (g.dat$start_position - pad)){
      cur.data <- data[data[,chr] == g.dat$chr & data[,pos] >= min.pos & data[,pos] <= max.pos, ]  
    }
    
    if (max.pos > (g.dat$end_position + pad)){
      cur.data <- data[data[,chr] == g.dat$chr & data[,pos] >= min.pos & data[,pos] <= max.pos, ]  
    }
    
    
    # make snp track
    snp.track <- DataTrack(
      data = -log10(cur.data[,pval]), 
      start = cur.data[,pos], 
      end = cur.data[,pos], 
      chromosome = cur.data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Trait -log10(p)", 
      background.title="orangered4",
      col = "orangered4",
      fill = "orangered4",
      col.frame="orangered4",
      frame = TRUE
    )
    
    if (nrow(cur.gtex) != 0){
      snp.gtex.track <- HighlightTrack(
        trackList = list(gtex.track, snp.track),
        start = cur.gtex$pos - 2, 
        end = cur.gtex$pos + 2,
        chromosome = cur.data[1,chr],
        col="lightgrey", 
        fill="lightgrey", 
        alpha=0.7
      )
    }
    
    
    if (gen == "hg38"){
      mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    } else {
      mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    }
    
    # gene track
    gtrack <- BiomartGeneRegionTrack(
      start = min.pos, 
      end = max.pos, 
      biomart = mart, 
      strand = "+-", 
      genome = gen, 
      chromosome = cur.data[1,chr], 
      name = "Genomic context", 
      showId = TRUE, 
      geneSymbol = FALSE,
      collapseTranscripts = TRUE,
      background.title = "cadetblue4"
    )
    
    # highlight the current gene
    hightlight.track <- HighlightTrack(
      trackList = list(gtrack),
      start = genes[genes$chr == cur.data[1,chr],]$start_position, 
      end = genes[genes$chr == cur.data[1,chr],]$end_position,
      chromosome = cur.data[1,chr]
    )
    
    # hmm track
    if (hmm.file != "NA"){
      # subset hmm
      cur.hmm <- hmm[hmm$chr == cur.data[1,chr] & hmm$end >= min.pos & hmm$start <= max.pos,]
      hmm.track <- AnnotationTrack(
        name = hmm.name, 
        start = cur.hmm$start, 
        end = cur.hmm$end, 
        id = cur.hmm$state,
        feature = cur.hmm$color,
        chromosome = cur.hmm$chr, 
        genome = gen,
        stacking="dense",
        col="transparent",
        groupAnnotation="feature",
        "#FFFFFF"="#FFFFFF",
        "#FF0000"="#FF0000",
        "#CD5C5C"="#CD5C5C",
        "#FFC34D"="#FFC34D",
        "#BDB76B"="#BDB76B",
        "#C0C0C0"="#C0C0C0",
        "#E9967A"="#E9967A",
        "#808080"="#808080",
        "#FFFF00"="#FFFF00",
        "#8A91D0"="#8A91D0",
        "#006400"="#006400",
        "#FF4500"="#FF4500",
        "#008000"="#008000",
        "#C2E105"="#C2E105",
        "#66CDAA"="#66CDAA",
        "#32CD32"="#32CD32"
      )
    }
    
    # scale track
    axTrack <- GenomeAxisTrack(genome = gen , chromosome = cur.data[1,chr])
    
    # idiogram track
    idxTrack <- IdeogramTrack(genome = gen , chromosome = cur.data[1,chr])
    
    # open png
    png(out.file, height = 8, width = 11, units = "in", res = 400)
    
    # plot it all
    if (hmm.file == "NA"){
      if (nrow(cur.gtex) == 0){
        plotTracks(
          c(idxTrack, axTrack, snp.track, hightlight.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 4),
          from = min.pos,
          to = max.pos
        )
      } else {
        plotTracks(
          c(idxTrack, axTrack, snp.gtex.track, hightlight.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 5, 4),
          from = min.pos,
          to = max.pos
        )
      }
      
    } else {
      if (nrow(cur.gtex) == 0){
        plotTracks(
          c(idxTrack, axTrack, snp.track, hightlight.track, hmm.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 4,1),
          from = min.pos,
          to = max.pos
        )
      } else {
        plotTracks(
          c(idxTrack, axTrack, snp.gtex.track, hightlight.track,hmm.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 5, 4,1),
          from = min.pos,
          to = max.pos
        )
      }
    }
    dev.off()
  }
}

plotTraitIndexGtex <- function(rs.file, data.file, gtex.file, pthresh = 1, chr = "chr", pos = "pos", pval = "P-value", rs.id.col = "SNP_id", gene.id.col = "GeneSymbol", gene.id.type = "hgnc", out.pref = "", gen = "hg19", pad = 500000, hmm.file = "NA", hmm.name = "ChromHMM"){
  
  # load my required library
  library(data.table)
  library(GenomicRanges)
  library(Gviz)
  library(biomaRt)
  
  # load genes
  rsids <- loadRSCoord(rs.file)
  
  # if hhm, load
  if (hmm.file != "NA"){
    hmm <- loadRoadmapHMM(hmm.file)  
  } else {
    hmm <- c()
  }
  
  # load gtex data
  gtex <- loadGTEx(gtex.file, gene.id.type)
  
  # load variant data
  data <- loadChargeRes(data.file, pval = pval, pt = pthresh)
  
  # counter
  i <- 1
  num.rsids <- length(unique(rsids[,gene.id.col]))
  
  # loop through and plot each desired group
  for (r in unique(rsids[,rs.id.col])){
    print(paste(as.character(i), as.character(num.rsids), sep = "/"))
    i <- i+1
    
    # subset to just this gene
    r.dat <- rsids[rsids[,rs.id.col] == r,][1,]
    
    
    min.pos <- r.dat$pos - pad
    min.pos <- r.dat$pos + pad
    
    
    # output file name
    out.file <- paste(out.pref, r, "png", sep = ".")  
    
    # subset to variants in group
    cur.data <- data[data[,chr] == r.dat$chr & data[,pos] >= min.pos & data[,pos] <= max.pos, ]
    
    if (nrow(cur.data) == 0) next
    
    # subset gtex to gene
    cur.gtex <- gtex[gtex$gene_id == r.dat[,gene.id.col],]
    
    if (nrow(cur.gtex) == 0){
      # y limit 
      ylim <- c(floor(min(-log10(cur.data[,pval]), na.rm = T)), ceiling(max(-log10(cur.data[,pval]), na.rm = T)))
    } else {
      # set out min and max positions to be updated as we go
      min.pos.new <- min(min.pos, cur.gtex$pos - 2000)
      max.pos.new <- max(max.pos, cur.gtex$pos + 2000)
      
      # y limit 
      ylim <- c(floor( min(min(-log10(cur.data[,pval]), na.rm = T),min(-log10(cur.gtex$pval_nominal), na.rm = T) ) ), ceiling( max(max(-log10(cur.data[,pval]), na.rm = T), max(-log10(cur.gtex$pval_nominal), na.rm = T))))
      
      # make snp track
      gtex.track <- DataTrack(
        data = -log10(cur.gtex$pval_nominal), 
        start = cur.gtex[,pos], 
        end = cur.gtex[,pos], 
        chromosome = cur.data[1,chr], 
        genome = gen, 
        ylim = ylim,
        name = paste(r.dat[,gene.id.col], "eQTL -log10(p)", sep = " "), 
        background.title="blue",
        col = "blue",
        pch = 17,
        col.frame="blue",
        frame = TRUE
      )
      
    }
    
    # fix boundaries for snp track
    if (min.pos < min.pos.new){
      min.pos <- min.pos.new
      cur.data <- data[data[,chr] == r.dat$chr & data[,pos] >= min.pos & data[,pos] <= max.pos, ]  
    }
    
    if (max.pos > max.pos.new){
      max.pos <- max.pos.new
      cur.data <- data[data[,chr] == r.dat$chr & data[,pos] >= min.pos & data[,pos] <= max.pos, ]  
    }
    
    
    # make snp track
    snp.track <- DataTrack(
      data = -log10(cur.data[,pval]), 
      start = cur.data[,pos], 
      end = cur.data[,pos], 
      chromosome = cur.data[1,chr], 
      genome = gen, 
      ylim = ylim,
      name = "Trait -log10(p)", 
      background.title="orangered4",
      col = "orangered4",
      fill = "orangered4",
      col.frame="orangered4",
      frame = TRUE
    )
    
    if (nrow(cur.gtex) != 0){
      snp.gtex.track <- HighlightTrack(
        trackList = list(gtex.track, snp.track),
        start = r.dat$pos - 2, 
        end = r.dat$pos + 2,
        chromosome = cur.data[1,chr],
        col="lightgrey", 
        fill="lightgrey", 
        alpha=0.7
      )
    } else {
      snp.highlight.track <- HighlightTrack(
        trackList = list(snp.track),
        start = r.dat$pos - 2, 
        end = r.dat$pos + 2,
        chromosome = cur.data[1,chr],
        col="lightgrey", 
        fill="lightgrey", 
        alpha=0.7
      )
    }
    
    
    if (gen == "hg38"){
      mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    } else {
      mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    }
    
    # gene track
    gtrack <- BiomartGeneRegionTrack(
      start = min.pos, 
      end = max.pos, 
      biomart = mart, 
      strand = "+-", 
      genome = gen, 
      chromosome = cur.data[1,chr], 
      name = "Genomic context", 
      showId = TRUE, 
      geneSymbol = FALSE,
      collapseTranscripts = TRUE,
      background.title = "cadetblue4"
    )
    
    # hmm track
    if (hmm.file != "NA"){
      # subset hmm
      cur.hmm <- hmm[hmm$chr == cur.data[1,chr] & hmm$end >= min.pos & hmm$start <= max.pos,]
      hmm.track <- AnnotationTrack(
        name = hmm.name, 
        start = cur.hmm$start, 
        end = cur.hmm$end, 
        id = cur.hmm$state,
        feature = cur.hmm$color,
        chromosome = cur.hmm$chr, 
        genome = gen,
        stacking="dense",
        col="transparent",
        groupAnnotation="feature",
        "#FFFFFF"="#FFFFFF",
        "#FF0000"="#FF0000",
        "#CD5C5C"="#CD5C5C",
        "#FFC34D"="#FFC34D",
        "#BDB76B"="#BDB76B",
        "#C0C0C0"="#C0C0C0",
        "#E9967A"="#E9967A",
        "#808080"="#808080",
        "#FFFF00"="#FFFF00",
        "#8A91D0"="#8A91D0",
        "#006400"="#006400",
        "#FF4500"="#FF4500",
        "#008000"="#008000",
        "#C2E105"="#C2E105",
        "#66CDAA"="#66CDAA",
        "#32CD32"="#32CD32"
      )
    }
    
    # scale track
    axTrack <- GenomeAxisTrack(genome = gen , chromosome = cur.data[1,chr])
    
    # idiogram track
    idxTrack <- IdeogramTrack(genome = gen , chromosome = cur.data[1,chr])
    
    # open png
    png(out.file, height = 8, width = 11, units = "in", res = 400)
    
    # plot it all
    if (hmm.file == "NA"){
      if (nrow(cur.gtex) == 0){
        plotTracks(
          c(idxTrack, axTrack, snp.highlight.track, gtrack),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 4),
          from = min.pos,
          to = max.pos
        )
      } else {
        plotTracks(
          c(idxTrack, axTrack, snp.gtex.track, gtrack),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 5, 4),
          from = min.pos,
          to = max.pos
        )
      }
      
    } else {
      if (nrow(cur.gtex) == 0){
        plotTracks(
          c(idxTrack, axTrack, snp.highlight.track, gtrack, hmm.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 4,1),
          from = min.pos,
          to = max.pos
        )
      } else {
        plotTracks(
          c(idxTrack, axTrack, snp.gtex.track, gtrack, hmm.track),
          showTitle = TRUE,
          sizes=c(1, 2, 5, 5, 4,1),
          from = min.pos,
          to = max.pos
        )
      }
    }
    dev.off()
  }
}