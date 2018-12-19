task subsetGds {
	Array[File] gds_files
	File variant_file
	String out_pref

	Int disk
	Int memory
	
	command {
		R --vanilla --args ${sep="," gds_files} ${variant_file} ${out_pref} < /variantResults/subsetGds.R
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

workflow w_subsetGds {
	# inputs
	Array[File] these_gds_files
	File this_variant_file
	String this_out_pref
	
	# other inputs
	Int this_memory
	Int this_disk

	
	call subsetGds {
		input: gds_files = these_gds_files, variant_file = this_variant_file, out_pref = this_out_pref, disk = this_disk, memory = this_memory
	}

	output {
		File subsetGds_file = subsetGds.out_file
	}
}