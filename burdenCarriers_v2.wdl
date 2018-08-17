task burdenCarriers {
	Array[File] gds_files
	File variant_file
	String out_pref
	File? sample_file

	Int disk
	Int memory
	
	command {
		R --vanilla --args ${sep="," gds_files} ${variant_file} ${out_pref} ${default="NA" sample_file} < /variantResults/burdenCarriers_v2.R
	}

	runtime {
		docker: "manninglab/variantresults:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File out_file = "${out_pref}.gds"
	}
}

workflow w_burdenCarriers {
	# inputs
	Array[File] these_gds_files
	File this_variant_file
	String this_out_pref
	File? this_sample_file

	# other inputs
	Int this_memory
	Int this_disk

	
	call burdenCarriers {
		input: gds_files = these_gds_files, variant_file = this_variant_file, out_pref = this_out_pref, sample_file = this_sample_file, disk = this_disk, memory = this_memory
	}

	output {
		File burdenCarriers_file = burdenCarriers.out_file
	}
}