# EPCombinations

Repository to accompany Martinez-Ara et al., 2022.

This repository contains all the scripts used to preprocess raw data to barcode counts, process barcode counts to activities and boost inidices and to analyse the processed data.

All scripts and pipelines were developed by Federico Comoglio and Miguel Martinez-Ara.

Structure:

* UpstreamPipeline contains a snakemake pipeline used to preprocess Upstream assay sequencing data.
* DownstreamPipeline contains a snakemake pipeline used to preprocess Downstream assay sequencing data.
* datteRo_1.2.1.tar.gz contains a custom R package developed by Federico Comoglio. The functions contained in this package are needed for processing of data in R scripts.
* Rproj contains all R scripts used to process intermediate count tables to activities and boost indices, to perform analyses and to generate figures.
* data contains processed activity and boost indices files, data needed for the design of cCRE-P libraries and external processed data needed for analyses (Except MicroC data and data from Bergman et al., 2021)


## Snakemake pipelines

Both pipelines contain 2 workflows for data processing. 
* ipcr worflow extracts barcodes and fragment identities associated to them from paired end iPCR sequencing data.
* ecn worflow extracts barcode counts from single end cDNA or pDNA sequencing data.

To run these pipelines config/run files have to be edited to point to the correct locations. The pipeline has to be run on the folder where Raw data is stored.

For iPCR data one config file and pipeline run is needed per locus design. A bowtie index is generated based on the design sequences stored in the data folder. A metadata file is needed for identification of samples, this is a tsv file containing a column "sample.id" that contains the name of the files to be processed until the read suffix ("_R1" or "_R2").

For cDNA/pDNA data one config file and pipeline run is needed per experiment/set of biological replicates.

See notes on example config files for example commands to run the pipelines.

## datteRo

To install datteRo use the following command in R.

````
install.packages("datteRo_1.2.1.tar.gz")
````

This package contains functions needed for the processing of data output from the snakemake pipelines.

## Rproj - analyses

Included in this folder:
 
* Postprocessing scripts from cDNA and pDNA counts to activities.
* Scripts to calculate (rescaled) boost indices for each of the datasets.
* Figures_and_analyses folder contains all scripts used for analyses shown in the paper and its figures.


## data

All processed data used for analyses is included here. In the folder 'Curated_Natoli' the updated TF motif database from Diaferia et al. (2016, 10.15252/embj.201592404) can be found. 'Housekeeping' contains the housekeeping classification of the promoters tested. 'tad_deconstr_desing' contains the TAD and DSH coordinates used to desing the libraries, as well as other external data used for analyses. This folder also includes the coordinates and sequences of the cCREs used for the libraries and their TF motif composition based on the TF motif database (fimo folders).

MicroC data used for analysis is not provided as it is too heavy for this repository. This can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130275. The data from Bergman et al., (2021) is available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184426. 

Intermediate files generated as output from the pipeline are not provided.




