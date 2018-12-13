# Functions to create paper-ready tables & figures of genome-wide sequence association studies

# variables that need to be defined:
#  google.bucket.loc <- "gs:/.../"

lapply(c("qqman","data.table","tools","RColorBrewer"), library, character.only = TRUE)

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
    print(head(assoc.data))

}