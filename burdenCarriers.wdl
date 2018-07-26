task burdenCarriers {
	File gds_file
	File assoc_file
	File group_file
	String out_pref
	Float? pval_thresh

	Int disk
	Int memory
	
	command {
		R --vanilla --args ${gds_file} ${assoc_file} ${group_file} ${out_pref} ${default="0.00001" pval_thresh} < /variantResults/burdenCarriers.R
	}

	runtime {
		docker: "manninglab/variantresults:0.1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		Array[File] carriers = glob("${out_pref}.*.carriers.csv")
		Array[File] genotypes = glob("${out_pref}.*.genotypes.csv")
	}
}

workflow w_burdenCarriers {
	# inputs
	Array[File] these_gds_files
	File this_assoc_file
	File this_group_file
	String this_out_pref
	Float? this_pval_thresh

	# other inputs
	Int this_memory
	Int this_disk

	scatter (this_gds_file in these_gds_files) {
		call burdenCarriers {
			input: gds_file = this_gds_file, assoc_file = this_assoc_file, group_file = this_group_file, out_pref = this_out_pref, pval_thresh = this_pval_thresh, disk = this_disk, memory = this_memory
		}
	}	
}