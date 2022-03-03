##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FC@NKI
# Design primers in batch
# (currently for DAT designs)
#
# For parametrization, see Primer3 manual (params explained): http://primer3.sourceforge.net/primer3_manual.htm#PRIMER_WT_TM_GT
# Python module for Primer3: https://libnano.github.io/primer3-py/quickstart.html#primer-design
#
# Background info:
#
# Output of primer3 contains header with following keys (n=7)
# PRIMER_LEFT_EXPLAIN, PRIMER_RIGHT_EXPLAIN, PRIMER_PAIR_EXPLAIN, PRIMER_LEFT_NUM_RETURNED, PRIMER_RIGHT_NUM_RETURNED, PRIMER_INTERNAL_NUM_RETURNED, PRIMER_PAIR_NUM_RETURNED
# (toggle/param PRIMER_EXPLAIN_FLAG appears not to work)
#
# SAM format specification:
# flag 99 = proper pair, first in pair, mate in reverse orientation (i.e. selects forward primer)
# flag 355 = as above, but not primary alignment
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load modules
import os
import sys
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import primer3 #primer3-py

# input fasta file name
# argv[0] is script file name
fa_fn = sys.argv[1]

# input bt2 index
bt2_index = sys.argv[2]

# output list
primer_sets = []

# output fasta file names (fwd, rev)
bn                 = os.path.basename(fa_fn)
out_fa_fwd         = os.path.splitext(bn)[0] + '_fwd.fasta'
out_fa_rev         = os.path.splitext(bn)[0] + '_rev.fasta'
out_fa_order_fwd   = os.path.splitext(bn)[0] + '_order_fwd.fasta'
out_fa_order_rev   = os.path.splitext(bn)[0] + '_order_rev.fasta'
out_sam            = os.path.splitext(bn)[0] + '_primer_aln.sam'
out_valid          = os.path.splitext(bn)[0] + '_primer_aln_valid.sam'
out_uniq           = os.path.splitext(bn)[0] + '_primer_aln_uniq.sam'
out_order          = os.path.splitext(bn)[0] + '_primer_aln_order.sam'

# iterate over input fasta
for record in SeqIO.parse(fa_fn, "fasta"):

    # generate input for primer3: ({id+seq},{params})
    specs      = {'SEQUENCE_ID': record.id,
                  'SEQUENCE_TEMPLATE': str(record.seq)}
    params     = {'PRIMER_NUM_RETURN': 10,
                  'PRIMER_OPT_SIZE': 27,
                  'PRIMER_PICK_INTERNAL_OLIGO': 0,
                  'PRIMER_MIN_SIZE': 20, #23
                  'PRIMER_MAX_SIZE': 35,
                  'PRIMER_OPT_TM': 65.0,
                  'PRIMER_MIN_TM': 62.0,
                  'PRIMER_MAX_TM': 68.0,
                  'PRIMER_PAIR_MAX_DIFF_TM': 3,
                  'PRIMER_GC_CLAMP': 1,
                  'PRIMER_MIN_GC': 20.0,
                  'PRIMER_OPT_GC_PERCENT': 50.0,
                  'PRIMER_MAX_GC': 80.0,
                  'PRIMER_MAX_POLY_X': 3,
                  'PRIMER_SALT_MONOVALENT': 10.0, # from MyTaq RedMix
                  'PRIMER_SALT_DIVALENT': 3.0,    # from MyTaq RedMix
                  'PRIMER_DNA_CONC': 400.0,       # 4ul of 1M in 10ul
                  'PRIMER_MAX_NS_ACCEPTED': 0,
                  'PRIMER_MAX_SELF_ANY': 10,
                  'PRIMER_MAX_SELF_END': 5,
                  'PRIMER_PAIR_MAX_COMPL_ANY': 10,
                  'PRIMER_PAIR_MAX_COMPL_END': 5,
                  'PRIMER_PAIR_WT_DIFF_TM': 1,
                  'PRIMER_TM_FORMULA': 1,
                  'PRIMER_PRODUCT_SIZE_RANGE': [[380,450]]}

    # design primer pairs
    this_set = primer3.bindings.designPrimers(specs, params)

    # add element id to primer pairs
    this_set['ID'] = record.id

    # if no primers found, output length is 8 (7 + id)
    if not len(this_set) == 8:
        primer_sets.append(this_set)

# number of elements with primers available
n_designed = str(len(primer_sets))
print('::Designed primers for ' + n_designed + ' input sequences')

# init output list
out_fwd = []
out_rev = []

# iterate over primer sets
for set in primer_sets:

    # compute number of pairs for given element
    n_pairs = int((len(set) - 8) / 22)

    # generate seq id, extract seq, append to output list
    for k in range(n_pairs):
        id_fwd = set['ID'] + '_' + str(k) + '_fwd'
        id_rev = set['ID'] + '_' + str(k) + '_rev'
        seq_fwd = Seq(set['PRIMER_LEFT_' + str(k) + '_SEQUENCE'])
        seq_rev = Seq(set['PRIMER_RIGHT_' + str(k) + '_SEQUENCE'])

        record_fwd = SeqRecord(seq_fwd, id_fwd, '', '')
        out_fwd.append(record_fwd)
        record_rev = SeqRecord(seq_rev, id_rev, '', '')
        out_rev.append(record_rev)

# write output fasta file
SeqIO.write(out_fwd, out_fa_fwd, 'fasta')
SeqIO.write(out_rev, out_fa_rev, 'fasta')

# PE primer alignment with bowtie2
cmd = 'nice -n 19 bowtie2 -f -t -p 6 -X 450 -k 2 --very-sensitive --no-unal --no-hd -x ' + bt2_index + ' -1 ' + out_fa_fwd + ' -2 ' + out_fa_rev + ' -S ' + out_sam
os.system(cmd)

# filter concordant pairs (retain valid primer sets)
cmd = "awk -F '\t' '{ if(($2 == 99) || ($2 == 355)) { print } }' " + out_sam + ' | sort > ' + out_valid
os.system(cmd)

# filter out primer sets with multiple aln
cmd = 'awk \'n=x[$1]{print n"\\n"$0;} {x[$1]=$0;}\' ' + out_valid + ' > tmp.txt; comm -23 ' + out_valid + ' tmp.txt | sort -V > ' + out_uniq + ' ; rm tmp.txt'
os.system(cmd)

# select best hit (i.e. lowest penalty, by id)
cmd = "cut -d'_' -f1 " + out_uniq + ' | paste - ' + out_uniq + ' | sort -u -k1,1 | sort -V > ' + out_order
os.system(cmd)

# select primer sequences to order
cmd = 'cut -f2 ' + out_order + " | sed 's/_fwd//g' > tmp.txt; grep -f tmp.txt -A 1 " + out_fa_fwd + " | grep -v '\-\-' > " + out_fa_order_fwd
os.system(cmd)
cmd = 'grep -f tmp.txt -A 1 ' + out_fa_rev + " | grep -v '\-\-' > " + out_fa_order_rev + '; rm tmp.txt'
os.system(cmd)
