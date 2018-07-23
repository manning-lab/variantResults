FROM r-base:3.5.0

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN echo "install.packages(c('data.table','dplyr','tidyr'), repos='http://cran.us.r-project.org')" > install.R && \
	echo "source('https://bioconductor.org/biocLite.R')" >> install.R && \
	echo "biocLite(c('SeqArray','SeqVarTools'))" >> install.R && \
	R --vanilla < install.R && \
	rm install.R

RUN cd / && \
	git clone https://github.com/manning-lab/variantResults.git && \
	cd variantResults && \
	git pull origin master