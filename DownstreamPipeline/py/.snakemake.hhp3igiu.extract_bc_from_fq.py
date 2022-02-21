
######## Snakemake header ########
import sys; sys.path.insert(0, "/DATA/usr/m.martinez.ara/miniconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05XR\x00\x00\x00/DATA/projects/epmatrix/mouse/DAT_EP010_EP013/ecn/Raw/pdna_E49_dat_br2_R1.fastq.gzq\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x05\x00\x00\x00readsq\nK\x00N\x86q\x0bsh\nh\x06ubX\x06\x00\x00\x00outputq\x0ccsnakemake.io\nOutputFiles\nq\r)\x81q\x0eX%\x00\x00\x00Results/pdna_E49_dat_br2_barcodes.tsvq\x0fa}q\x10(h\x08}q\x11X\x03\x00\x00\x00tsvq\x12K\x00N\x86q\x13sh\x12h\x0fubX\x06\x00\x00\x00paramsq\x14csnakemake.io\nParams\nq\x15)\x81q\x16}q\x17h\x08}q\x18sbX\t\x00\x00\x00wildcardsq\x19csnakemake.io\nWildcards\nq\x1a)\x81q\x1bX\x10\x00\x00\x00pdna_E49_dat_br2q\x1ca}q\x1d(h\x08}q\x1eX\x03\x00\x00\x00ecnq\x1fK\x00N\x86q sX\x03\x00\x00\x00ecnq!h\x1cubX\x07\x00\x00\x00threadsq"K\x01X\t\x00\x00\x00resourcesq#csnakemake.io\nResources\nq$)\x81q%(K\x01K\x01e}q&(h\x08}q\'(X\x06\x00\x00\x00_coresq(K\x00N\x86q)X\x06\x00\x00\x00_nodesq*K\x01N\x86q+uh(K\x01h*K\x01ubX\x03\x00\x00\x00logq,csnakemake.io\nLog\nq-)\x81q.X5\x00\x00\x00logs/py/pdna_E49_dat_br2_extract_barcodes_from_fq.logq/a}q0h\x08}q1sbX\x06\x00\x00\x00configq2}q3(X\x0b\x00\x00\x00descriptionq4X"\x00\x00\x00EP-SuRE ECN data analysis pipelineq5X\x07\x00\x00\x00datasetq6X%\x00\x00\x00ECN runs (pDNA + cDNA), SE sequencingq7X\x08\x00\x00\x00metadataq8ccollections\nOrderedDict\nq9)Rq:(X\x0b\x00\x00\x00mate_suffixq;X\x03\x00\x00\x00_R1q<X\x0c\x00\x00\x00input_formatq=X\x08\x00\x00\x00fastq.gzq>uX\x03\x00\x00\x00dirq?h9)Rq@(X\t\x00\x00\x00fastq_dirqAX5\x00\x00\x00/DATA/projects/epmatrix/mouse/DAT_EP010_EP013/ecn/RawqBX\x07\x00\x00\x00resultsqCX\x07\x00\x00\x00ResultsqDX\x04\x00\x00\x00logsqEX\x04\x00\x00\x00LogsqFX\t\x00\x00\x00conda_envqGXC\x00\x00\x00/home/m.martinez.ara/mydata/GitLab/epmatrix/epsure_pipeline/v2/envsqHX\n\x00\x00\x00snakerulesqIXD\x00\x00\x00/home/m.martinez.ara/mydata/GitLab/epmatrix/epsure_pipeline/v2/rulesqJuuX\x04\x00\x00\x00ruleqKX\x18\x00\x00\x00extract_barcodes_from_fqqLub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FC@NKI
# EP-SuRE pipeline

# Extract barcodes from fastq files (cDNA and pDNA SE data)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load modules
import pysam
import regex

# define seq immediately downstream of barcode (truncated to 30nt)
# allow for 2 mismatches (fuzzy search)
# note full seq: 'CCTAGCTAACTATAACGGTCCTAAGGTAGCGAAGGATCCATGCCC'
# (partly common to iPCR linker)
downstream_seq = '(' + 'CCTAGCTAACTATAACGGTCCTAAGGTAGC' + '){e<3}'

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
