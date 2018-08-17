# __variantResults__
Tools for exploring single and aggregate variant results. All workflows are written with a *Workflow definition language* (*WDL*) wrapper for transportability across compute environments.

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). Contributing authors include:

* Tim Majarian
* Alisa Manning

### General Dependencies

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Cromwell](http://cromwell.readthedocs.io/en/develop/)

# Workflows

## burdenCarriers

### Description 

This workflow creates a new, subsetted genotype file (in GDS format) for a defined set of variants. This new genotype file contains only those variants specified and only those samples that carry at least one minor allele in at least one variant. The initial motivation for this workflow is in exploring rare and low-frequency variant associations.

### R packages

All code and necessary packages are available through a [Docker image](https://hub.docker.com/r/manninglab/variantresults/) as well as through the [Github repository](https://github.com/manning-lab/variantResults/).

* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [dplyr](https://dplyr.tidyverse.org/)
* [tidyr](https://tidyr.tidyverse.org/)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

### Main Functions

#### burdenCarriers.R
This script generates the subsetted GDS file for the variants specified.

Inputs:
* gds_files : comma separated list of gds file paths (file, GDS)
* variant_file : file containing the variants to subset by with a minimum header of: chromosome, position, ref, alt (file, CSV/TSV)
	* the optional column *minor.allele* should specify either *ref* or *alt* for each variant, determining which allele should be treated as the minor. If no *minor.allele* column is specified, *alt* is assumed to be the minor.
* out_pref : prefix for output filename (string)
* *sample_file* : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (TXT)
* disk : amount of disk space in GB for task (INT)
* memory : amount of memory in GB for task (INT)

Outputs :
* out_file : genotype file containing only those variants specified in *variant_file* and only those samples who carrier a minor allele (GDS, *out_pref*.gds)

### Workflow Outputs

The only output is *out_file* from above
