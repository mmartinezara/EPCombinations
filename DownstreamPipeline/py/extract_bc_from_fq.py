##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MMA and FC@NKI
# P-E-SuRE pipeline (Split design)

# Extract barcodes from fastq files (cDNA and pDNA SE data)
# This script is based on the equivalent one of the EP-SuRE pipeline V2. The only difference is the sequence downstream of the barcode which is specific for each setup.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load modules
import pysam
import regex

# define seq immediately downstream of barcode (truncated to 30nt)
# allow for 2 mismatches (fuzzy search)
# note full seq: 'GCTAGGTTCGCTACCTTAGGACCGTTATAGTTATCATCGAGCTCG'
# (partly common to iPCR linker)
downstream_seq = '(' + 'GCTAGGTTCGCTACCTTAGGACCGTTATAG' + '){e<3}'

# open output file
tsv_out = open(snakemake.output["tsv"], "w")

# open input fastq stream
with pysam.FastxFile(snakemake.input["reads"]) as fq_in:

    # iterate over reads
    for read in fq_in:

        # extract read sequence
        seq = read.sequence

        # identify downstream seq position
        match = regex.search(downstream_seq, seq, regex.BESTMATCH)

        # if no match, skip to next
        if match is None:
            continue

        # extract barcode
        end_bc = match.span()[0]
        barcode = seq[0:end_bc]

        # if barcode intact and no N in barcode, write to file
        if((len(barcode) == 20) and ('N' not in barcode)):
            # write to output file
            tsv_out.write(barcode + '\n')

tsv_out.close()
