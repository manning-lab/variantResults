# subsetting single variant results by a bed file

# Inputs:
# this_results_file: csv of tsv with aggregate or single variants results, must have column headers including "chr" (or "chromosome") and "pos" (or "position")
# this_bed_file: csv or tsv in (at least) 3 column BED format with no column headers (columns in order of chromosome, start coordinate, end coordinate)
# this_out_pref: output file prefix (optional, default = "subset.results")
# this_memory: amount of memory in GB for running the task
# this_disk: amount of disk space in GB for running the task

# Outputs:
# subset_results: csv with only those results intersecting the intervals in this_bed_file, with the same columns as this_results_file (default = "subset.results.csv")
task subset {
	File results_file
	File bed_file
	String? out_pref
	String out_name = select_first([out_pref,"subset.results"])

	Int memory
	Int disk

	command <<<
		echo "lapply(c('data.table','GenomicRanges'),library, character.only = T)" > script.R && \
		echo "res.file <- '${results_file}'" >> script.R && \
		echo "bed.file <- '${bed_file}'" >> script.R && \
		echo "res.data <- fread(res.file, data.table = F, stringsAsFactors = F)" >> script.R && \
		echo "ntokeep <- names(res.data)" >> script.R && \
		echo "chrpos <- c('chr','pos')" >> script.R && \
		echo "if('chromosome' %in% ntokeep){chrpos[1] <- 'chromosome'}" >> script.R && \
		echo "if('position' %in% ntokeep){chrpos[2] <- 'position'}" >> script.R && \
		echo "res.data[,chrpos[1]] <- as.character(res.data[,chrpos[1]])" >> script.R && \
		echo "bed.data <- fread(bed.file, data.table = F, stringsAsFactors = F, header = F)" >> script.R && \
		echo "bed.data[,'V1'] <- as.character(bed.data[,'V1'])" >> script.R && \
		echo "if(startsWith(res.data[1,chrpos[1]],'chr') & !(startsWith(bed.data[1,'V1'], 'chr'))){bed.data[,'V1'] <- sub('^','chr',bed.data[,'V1'])}" >> script.R && \
		echo "if(!(startsWith(res.data[1,chrpos[1]],'chr')) & startsWith(bed.data[1,'V1'], 'chr')){bed.data[,'V1'] <- sub('chr','',bed.data[,'V1'])}" >> script.R && \
		echo "res.data[,'ident'] <- seq(1,nrow(res.data))" >> script.R && \
		echo "res.gr <- GRanges(seqnames = res.data[,chrpos[1]], ranges = IRanges(start = res.data[,chrpos[2]], end = res.data[,chrpos[2]]), ident = res.data[,'ident'])" >> script.R && \
		echo "bed.gr <- GRanges(seqnames = bed.data[,'V1'], ranges = IRanges(start = bed.data[,'V2'], end = bed.data[,'V3']))" >> script.R && \
		echo "res.ovp <- res.gr[queryHits(findOverlaps(res.gr,bed.gr)),]" >> script.R && \
		echo "res.data <- res.data[res.data[,'ident'] %in% as.data.frame(res.ovp)[,'ident'],ntokeep]" >> script.R && \
		echo "fwrite(res.data, file = paste0('${out_name}','.csv'), row.names = F, sep = ',')" >> script.R && \
		R --vanilla < script.R
	>>>

	runtime {
		docker: "manninglab/variantresults:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File out_file = "${out_name}.csv"
	}
}

workflow subsetResults {
	File this_results_file
	File this_bed_file
	String? this_out_pref
	Int this_memory
	Int this_disk

	call subset {
		input: results_file = this_results_file, bed_file = this_bed_file, out_pref = this_out_pref, memory = this_memory, disk = this_disk
	}

	output {
		File subset_results = subset.out_file
	}
}