FROM robbyjo/r-mkl-bioconductor:3.4.3-16.04-2018.1

RUN apt-get update && apt-get install -y git build-essential libcurl4-openssl-dev

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN echo "install.packages(c('data.table','dplyr','tidyr'), repos='http://cran.us.r-project.org')" > install.R && \
	echo "source('https://bioconductor.org/biocLite.R')" >> install.R && \
	echo "biocLite(c('SeqArray','SeqVarTools'))" >> install.R && \
	R --vanilla < install.R && \
	rm install.R

RUN cd / && \
	git clone https://github.com/manning-lab/variantResults.git && \
	cd variantResults && \
	git pull origin master && \
	echo thisHIHi