##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MMA & FC@NKI
# EP-SuRE pipeline v3

# Identify fragments, assign orientation and barcodes to iPCR
# We include a step to identify empty fragments so that we can use them to normalise across libraries. This is based in the same step introduced in the Split design pipeline.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load modules
import pysam
import regex
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# open bam file and output file
bam_in = pysam.AlignmentFile(snakemake.input["bam"], "rb")
tsv_out = open(snakemake.output["tsv"], "w")

# define mate 1 linker sequence + allow for 2 mismatches (fuzzy search)
linker_seq = '(' + 'CCTAGCTAACTATAACGGTCCTAAGGTAGCGAACCAGTGAT' + '){e<3}'

# define mate 2 linker sequence before barcode + allow for 2 mismatches (fuzzy search)
linker_seq2 = '(' + 'GACCGTTATAGTTAGCTAGG' + '){e<3}'
# define mate 2 linker sequence after barcode + allow for 2 mismatches (fuzzy search)
linker_seq3 = '(' + 'AGATCGGAAGAGCGTCG' + '){e<3}'

linker_seq4 = '(' + 'CTGGAGATCGGA' + '){e<3}'
# iterate over reads
for read in bam_in.fetch(until_eof=True):
     #identify empty combinations
  if read.is_paired and read.is_read2 and read.is_unmapped and (read.mate_is_unmapped):
    
    # read1 sequence
            seq = read.query_sequence
            
    # identify linker2 position before barcode read 2
            match = regex.search(linker_seq2, seq, regex.BESTMATCH)
            
    # identify linker3 position after barcode read 2

            match2 = regex.search(linker_seq3, seq, regex.BESTMATCH)
             # if no linker match, skip to next
            if match2 is not None and match is not None:
                  # extract fragments identity
                  frag1_chr = 'empty'
                  frag2_chr = 'empty'
                  frag1_strand = '*'
                  frag2_strand = '*'

                  # extract barcode
                  start_bc = match.span()[1]
                  end_bc = match2.span()[0]
                  barcode = Seq(seq[start_bc:end_bc], generic_dna)

                  # get reverse reverse_complement barcode

                  # if no N in barcode, write to file
                  if('N' not in barcode and len(barcode) > 0):
                    # write to output file
                    tsv_out.write(str(barcode.reverse_complement()) + '\t' + frag1_chr + '\t' + frag1_strand + '\t' + frag2_chr + '\t' + frag2_strand + '\n')
  else:
    if read.is_paired and read.is_read1 and read.is_unmapped and (read.mate_is_unmapped):
      # read1 sequence
      seq = read.query_sequence
      match3 = regex.search(linker_seq4, seq, regex.BESTMATCH)
      match = regex.search(linker_seq, seq, regex.BESTMATCH)
      if match3 is None and match is not None:
        # extract fragments identity
        frag1_chr = 'unknown'
        frag2_chr = 'unknown'
        frag1_strand = '*'
        frag2_strand = '*'
        #extract barcode
        end_bc = match.span()[0]
        barcode = seq[0:end_bc]
        # if no N in barcode, write to file
        if('N' not in barcode and len(barcode) > 0):
          # write to output file
          tsv_out.write(barcode + '\t' + frag1_chr + '\t' + frag1_strand + '\t' + frag2_chr + '\t' + frag2_strand + '\n')
  
    else:
      # require paired (100%) + mate1 + mate2 mapped + expected barcode+linker seq length
      # note: (4, 61) encodes 61S in CIGAR
      if read.is_paired and (not read.is_unmapped) and read.is_read1 and (not read.mate_is_unmapped) and (read.cigartuples[0] == (4, 61) or read.cigartuples[-1] == (4, 61)):
        # check if mate1 maps to reverse strand
        # note: if FLAG contains bit value 4 --> SEQ in SAM file is reverse complement
        if(read.is_reverse):

            # fragment 1 maps to '-'
            frag1_strand = '-'

            # reverse complement read sequence
            seq = Seq(read.query_sequence)
            seq = seq.reverse_complement()
            seq = str(seq)

            # identify linker position
            match = regex.search(linker_seq, seq, regex.BESTMATCH)

            # if no linker match, skip to next
            if match is None:
                continue

        else:
            # fragment 1 maps to '+'
            frag1_strand = '+'

            # read sequence is forward
            seq = read.query_sequence

            # identify linker position
            match = regex.search(linker_seq, seq, regex.BESTMATCH)

            # if no linker match, skip to next
            if match is None:
                continue

        # extract fragments identity
        frag1_chr = read.reference_name
        frag2_chr = read.next_reference_name

        # extract fragment 2 orientation
        if(read.mate_is_reverse):
            frag2_strand = '+'
        else:
            frag2_strand = '-'

        # extract barcode
        end_bc = match.span()[0]
        barcode = seq[0:end_bc]

        # if no N in barcode, write to file
        if('N' not in barcode and len(barcode) > 0):
            # write to output file
            tsv_out.write(barcode + '\t' + frag1_chr + '\t' + frag1_strand + '\t' + frag2_chr + '\t' + frag2_strand + '\n')

bam_in.close()
tsv_out.close()
