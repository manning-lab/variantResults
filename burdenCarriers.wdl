task burdenCarriers {
	File gds_file
	File variant_file
	String out_pref
	File? sample_file

	Int disk
	Int memory
	
	command {
		R --vanilla --args ${gds_file} ${variant_file} ${out_pref} ${default="NA" sample_file} < /variantResults/burdenCarriers.R
	}

	runtime {
		docker: "manninglab/variantresults:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File? out_file = select_first(glob("${out_pref}.*.gds"))
	}
}

task combineGds {
	Array[File?] gds_files
	String out_pref

	Int disk
	Int memory

	command {
		R --vanilla --args ${sep="," gds_files} ${out_pref} < /variantResults/combineGds.R
	}

	runtime {
		docker: "manninglab/variantresults:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File? full_gds = select_first(glob("${out_pref}.gds"))
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

	scatter (this_gds_file in these_gds_files) {
		call burdenCarriers {
			input: gds_file = this_gds_file, variant_file = this_variant_file, out_pref = this_out_pref, sample_file = this_sample_file, disk = this_disk, memory = this_memory
		}
	}

	call combineGds {
		input: gds_files = burdenCarriers.out_file, out_pref = this_out_pref, disk = this_disk, memory = this_memory
	}

	output {
		File? final_gds = combineGds.full_gds
	}
}