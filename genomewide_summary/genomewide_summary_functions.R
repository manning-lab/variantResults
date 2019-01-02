# Functions to create paper-ready tables & figures of genome-wide sequence association studies

get.data <- function(fname) {
   if(!file.exists(fname)) {
	     print(paste("Copying",saige.file,"from google bucket"))
	system(paste("gsutil cp ",google.bucket.loc,fname," ",fname,sep=""))
   }
   if(grepl(".gz$",fname)) {
            print("Reading in GZipped file")
            assoc.data <- fread(cmd=paste("gunzip -c ",fname,sep=""),data.table=T,key="chr,pos,ref,alt")
   } else {
            print("Reading in file")
            assoc.data <- fread(fname,data.table=T,key="chr,pos,ref,alt")
   }
}

# variables that need to be defined:
#  google.bucket.loc <- "gs:/.../"
#lapply(c("qqman","data.table","tools","RColorBrewer"), library, character.only = TRUE)

do.saige.data <- function(saige.file) {
    print(saige.file)
    if(!file.exists(saige.file)) {
        print(paste("Copying",saige.file,"from google bucket"))
        system(paste("gsutil cp ",google.bucket.loc,saige.file," ",saige.file,sep=""))
    }
    assoc.data <- fread(cmd=paste("gunzip -c",saige.file),data.table=T,showProgress=T)
    print(dim(assoc.data))
    #print(head(assoc.data))
    #print(tail(assoc.data))
    
    label <- saige.file

    assoc.data$MarkerName <- paste(assoc.data$CHR,"-",assoc.data$POS,"-",assoc.data$Allele1,"-",assoc.data$Allele2,sep="")
    
    assoc.data$minor.allele <- ifelse(assoc.data$AF_Allele2<0.5,"alt","ref")
    assoc.data$maf <- ifelse(assoc.data$AF_Allele2<0.5,assoc.data$AF_Allele2,1-assoc.data$AF_Allele2)
    assoc.data$mac <- ifelse(assoc.data$AF_Allele2<0.5,assoc.data$AC_Allele2,2*assoc.data$N-assoc.data$AC_Allele2)
    tmp <- sub("chr","",assoc.data$CHR)
    assoc.data$CHR <- as.numeric(ifelse(tmp=="X",23,tmp))
    
    setnames(assoc.data, c("CHR","POS","Allele1","Allele2","N","p.value","AF_Allele2"), 
         c("chr","pos","ref","alt","n","pvalue","ALTFreq"))    
    assoc.data <- assoc.data[,list(MarkerName,chr,pos,ref,alt,minor.allele,maf,mac,n,pvalue,SNPID,BETA,SE,ALTFreq)]

    print(table(assoc.data$chr))
    head(assoc.data)
    fwrite(assoc.data[pvalue<5e-6,],file=paste(label,"Plt5e-6.csv",sep="_"),quote=FALSE)
    fwrite(assoc.data[pvalue<0.01,],file=paste(label,"Plt0.01.csv",sep="_"),quote=FALSE)
    fwrite(assoc.data,file=paste(label,"allvariants.csv",sep="_"),quote=FALSE)

    system(paste("gsutil cp ",label,"_Plt5e-6.csv ",google.bucket.loc,sep=""))
    system(paste("gsutil cp ",label,"_Plt0.01.csv ",google.bucket.loc,sep=""))
    system(paste("gzip -f ",label,"_allvariants.csv" ,sep=""))
    system(paste("gsutil cp ",label,"_allvariants.csv.gz ",google.bucket.loc,sep=""))

	fwrite(assoc.data[pvalue<0.01,],file=paste(label,"ZOOM.txt",sep="_"),quote=FALSE,sep="\t")
	system(paste("htslib-1.9/bgzip -f ",label,"_ZOOM.txt",sep=""))
	system(paste("htslib-1.9/tabix -f -b 3 -e 3 -s 2 --skip 1 ",label,"_ZOOM.txt.gz",sep=""))
	system(paste("gsutil cp ",label,"_ZOOM.txt.gz ",label,"_ZOOM.txt.gz.tbi ",google.bucket.loc,sep=""))

}

### Create BINS (500 MB windows)


do.big.summary<- function(i,assoc.data.index.results){
    #rbind(
    tmp1 <- cbind(assoc.data.index.results[[i]][["index.result"]],
                  assoc.data.index.results[[i]][["assoc.data.index"]][,.(minpos=min(pos),maxpos=max(pos),nvars=length(pos))])
    tmp3 <- cbind(tmp1,assoc.data.index.results[[i]][["assoc.data.index"]])
    return(tmp3)
    #      assoc.data.pooled.index.results[[i]][["assoc.data.index"]])
}
get.binned <- function(assoc.data) {
    assoc.data <- assoc.data[order(chr,pos)]
    print(head(assoc.data))
    print(paste("Nrows remaining",nrow(assoc.data)))

    list.results <- list()
    
    i<-1
    
    while(nrow(assoc.data)>0) {
        print(paste("Index",i))
        list.results[[i]] <- list()
        list.results[[i]][["index.result"]] <- assoc.data[which.min(pvalue),]

        index.chr <- unlist(list.results[[i]][["index.result"]][1,"chr"])
        print(index.chr)
        index.pos <- unlist(list.results[[i]][["index.result"]][1,"pos"])
        print(index.pos)
        list.results[[i]][["assoc.data.index"]] <- assoc.data[chr == index.chr & (pos > (index.pos - 500000) &  pos < (index.pos+500000))]
        i <- i+1                                                 
        assoc.data <- assoc.data[chr != index.chr | (chr == index.chr) & (pos < (index.pos - 500000) |  pos > (index.pos+500000))]
                                
        print(head(assoc.data))
        print(paste("Nrows remaining",nrow(assoc.data)))
    }
    return(list.results)
}

do.it.all <- function(in.file,out.file) {
    assoc.data <- get.data(in.file)
    assoc.data <- assoc.data[mac>20]
    assoc.data$OR =  exp(assoc.data$BETA)

    head(assoc.data)
    summary(assoc.data)
    assoc.data.index.results <- get.binned(assoc.data)

    assoc.data.index.results.summary <- t(sapply(1:length(assoc.data.index.results),
                                                            function(i){cbind(assoc.data.index.results[[i]][["index.result"]],
                                                                              assoc.data.index.results[[i]][["assoc.data.index"]][,.(minpos=min(pos),maxpos=max(pos),nvars=length(pos))])}))
    assoc.data.index.results.summary

    assoc.data.index.results.summary.all <- rbindlist(sapply(1:length(assoc.data.index.results),
                                                                do.big.summary,
                                                                assoc.data.index.results=assoc.data.index.results,
                                                                simplify = FALSE, USE.NAMES = FALSE))

    fwrite(data.frame(assoc.data.index.results.summary),file=paste(out.file,".binned.csv",sep=""))
    fwrite(data.frame(assoc.data.index.results.summary.all),file=paste(out.file,".binned.all.csv",sep=""))
    system(paste("gsutil cp ",out.file,".binned.all.csv ",out.file,".binned.csv ",google.bucket.loc,sep=""))
}


### Manhattan Plots

get.ancestry <- function(anc.label) {

    assoc.data <- fread(get.data(assoc.names[[anc.label]]),data.table=T,showProgress=T)
    print(dim(assoc.data))
 
    print(table(p.lt.0.01=assoc.data$pvalue<0.01))
    print(table(p.lt.0.01=assoc.data$pvalue<0.01,common=assoc.data$maf>=0.01))
    print(table(p.lt.0.01=assoc.data$pvalue<0.01,rare=assoc.data$maf<0.01 & assoc.data$mac>20))
    print(table(p.lt.0.01=assoc.data$pvalue<0.01,ultrarare=assoc.data$mac<=20))
        
    print(paste("Min p-value:",min(assoc.data$pvalue)))

  
	return(assoc.data)
}

do.ancestry.manhattan <- function(anc.label,assoc.data.all,maxy=10) {
  	options(repr.plot.width=12, repr.plot.height=4)

    print(paste("Common SNPs:",sum(assoc.data.all$maf>=0.01)))
    #assoc.ps[["maf"]] >= 0.01]))
    do.manhattan(data.frame(assoc.data.all[which(assoc.data.all$maf>=0.01),c("chr","pos","pvalue")]),p="pvalue",ylim=c(2,maxy))
   
    print(paste("Rare SNPs:",sum(assoc.data.all$maf<0.01 & assoc.data.all$mac>20))) 
    #assoc.ps[["maf"]] < 0.01 & assoc.ps[["mac"]]>20])
    do.manhattan(data.frame(assoc.data.all[which(assoc.data.all$maf<0.01 & assoc.data.all$mac>20),c("chr","pos","pvalue")]),p="pvalue",ylim=c(2,maxy))

    print(paste("Ultrarare SNPs:",sum(assoc.data.all$mac<=20)))
    #assoc.ps[["mac"]]<=20])
    do.manhattan(data.frame(assoc.data.all[which(assoc.data.all$mac<=20),c("chr","pos","pvalue")]),p="pvalue",ylim=c(2,maxy))
}


do.manhattan <- function(assoc.data,ylim=c(3,32),chr="chr",bp="pos",p="P.value") {
    manhattan(assoc.data,chr=chr,bp=bp,p=p,
         suggestiveline = -log10(5e-07), genomewideline = -log10(5e-08) ,ylim=ylim)
}


do.manhattan.saige <- function(in.file,maxy=10) {
    assoc.data <- get.data(in.file)
    
    options(repr.plot.width=12, repr.plot.height=4)

    print(paste("Common SNPs:",sum(assoc.data$maf>=0.01)))
    do.manhattan(data.frame(assoc.data[which(assoc.data$maf>=0.01),
                                              c("chr","pos","pvalue")]),chr="chr",p="pvalue",ylim=c(2,maxy))
   
    print(paste("Rare SNPs:",sum(assoc.data$maf<0.01))) 
    do.manhattan(data.frame(assoc.data[which(assoc.data$maf<0.01),
                                              c("chr","pos","pvalue")]),chr="chr",p="pvalue",ylim=c(2,maxy))
}

## QQ plots
lam.new <- function(x,p=.5){
  x = x[!is.na(x)]
  #chisq <- qchisq(1-x,1)
  x.quantile <- quantile(x,p)
  round((qchisq(1-x.quantile,1)/qchisq(p,1)),2)
}

qqplot.saige <- function(in.file) {
    assoc.data <- get.data(in.file)
    head(assoc.data)
    assoc.ps <- list()

    assoc.ps[["ps"]] <- unlist(assoc.data$pvalue)
    assoc.ps[["maf"]] <- unlist(assoc.data$maf)
    rm(assoc.data)

    options(repr.plot.width=4, repr.plot.height=4)

    print(paste("Common SNPs:",sum(assoc.ps[["maf"]] >= 0.01)))
    print(summary(assoc.ps[["ps"]][assoc.ps[["maf"]] >= 0.01]))
    print(table(assoc.ps[["maf"]] >= 0.01))
    lam.meta.common <- lam.new(assoc.ps[["ps"]][assoc.ps[["maf"]] >= 0.01])
    print(lam.meta.common)
    qqpval3(assoc.ps[["ps"]][assoc.ps[["maf"]] >= 0.01],
        main=paste0('GC (MAF >= 1%)=',lam.meta.common))

    print(paste("Rare SNPs:",sum(assoc.ps[["maf"]] < 0.01 ) )) 
    summary(assoc.ps[["ps"]][assoc.ps[["maf"]] < 0.01 ])
    print(table(assoc.ps[["maf"]] < 0.01 ))
    lam.meta.rare <- lam.new(assoc.ps[["ps"]][assoc.ps[["maf"]] < 0.01 ])
    print(lam.meta.rare)
    qqpval3(assoc.ps[["ps"]][assoc.ps[["maf"]] < 0.01 ],
        main=paste0('GC (MAF < 1%)=',lam.meta.rare))
    rm(assoc.ps)
    gc()
}
