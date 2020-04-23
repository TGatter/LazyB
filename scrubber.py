#!/usr/bin/env python

import networkx as nx
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

paf_file = sys.argv[1] # input needs to be dual=yes !!!!
paf_illu_file = sys.argv[2]
seq_file_base = sys.argv[3]
output_file = sys.argv[4]
output = open(output_file, "w")

if seq_file_base.endswith("fa") or seq_file_base.endswith("fasta"):
    seq_dict = SeqIO.index_db(seq_file_base+".idx", seq_file_base, "fasta")
else:
    seq_dict = SeqIO.index_db(seq_file_base+".idx", seq_file_base, "fastq")

def find_sequence(seqid, start, end):
   return seq_dict[seqid].seq[start:end + 1]

def reverse_complement(seq):
   return seq.reverse_complement()

def write_sequence(f, seq):
    for i in range(0, len(seq), 60):
        f.write(str(seq[i:i + 60]) + "\n")


# Collect overlaps for all 
illu_to_ranges = {}

nano_set = set()

for line in open(paf_illu_file, "r"):

    linetemp = line.rstrip().split("\t")
    if len(linetemp) == 1:
        continue;

#0	qseqid	string 	Query sequence name         Illumina                 
#1	qlen 	int 	Query sequence length
#2	qstart	int 	Query start (0-based)
#3	qend	int 	Query end (1-based)
#4 	strand	char 	Relative strand: "+" or "-"

#5	sseqid 	string 	Target sequence name        Nano
#6	slen 	int 	Target sequence length
#7	sstart	int 	Target start on original strand (0-based)
#8	send	int 	Target end on original strand (0-based)
#9 	int 	Number of residue matches
#10 	int 	Alignment block length
#11 	int 	Mapping quality (0-255; 255 for missing) TODO

    id_1, id_2 = linetemp[0], linetemp[5]  # remember about 0-based counting

    len_1, len_2 = int(linetemp[1]), int(linetemp[6])  # remember about 0-based counting

    s_1, s_2 = int(linetemp[2]), int(linetemp[7])
    e_1, e_2 = int(linetemp[3]), int(linetemp[8])

    if e_1 - s_1 < 500:
        continue

    if id_2 not in illu_to_ranges:
        illu_to_ranges[id_2] = {}

    if id_1 not in illu_to_ranges[id_2]:
        illu_to_ranges[id_2][id_1] = (s_2, e_2)

    nano_set.add(id_2)

def print_corrected_read(cid, seq_to_ranges, illu_to_ranges, length):
    join = []
    for sid in seq_to_ranges:
        s, e, _ = seq_to_ranges[sid]
        join.append( (s,e) )
    for sid in illu_to_ranges:
        s, e = illu_to_ranges[sid]
        join.append( (s,e) )

    join.sort()

    covered = []
    for s, e in join:
        if len(covered) == 0:
            covered.append( (s, e) )
            continue

        cs = covered[-1][0]
        ce = covered[-1][1] 
        if cs <= e and s <= ce:
            covered[-1] = (min(s, cs), max(e, ce)) # join
        else:
            covered.append( (s, e) ) # add

    for i, (cs, ce) in enumerate(covered):
        #print cid, cs, ce
        output.write( ">" + cid + "_" + str(i) + "\n")
        write_sequence(output, find_sequence(cid, max(cs, 200), min(ce, length - 200)) )


# Collect overlaps for all 
seq_to_ranges = {}
cid = ""
length = 0

for line in open(paf_file, "r"):

    linetemp = line.rstrip().split("\t")
    if len(linetemp) == 1:
        continue;

#0	qseqid	string 	Query sequence name                          
#1	qlen 	int 	Query sequence length
#2	qstart	int 	Query start (0-based)
#3	qend	int 	Query end (1-based)
#4 	strand	char 	Relative strand: "+" or "-"

#5	sseqid 	string 	Target sequence name
#6	slen 	int 	Target sequence length
#7	sstart	int 	Target start on original strand (0-based)
#8	send	int 	Target end on original strand (0-based)
#9 	int 	Number of residue matches
#10 	int 	Alignment block length
#11 	int 	Mapping quality (0-255; 255 for missing) TODO

    id_1, id_2 = linetemp[0], linetemp[5]  # remember about 0-based counting
    if id_1 == id_2:
        continue
  
    len_1, len_2 = int(linetemp[1]), int(linetemp[6])  # remember about 0-based counting

    s_1, s_2 = int(linetemp[2]), int(linetemp[7])
    e_1, e_2 = int(linetemp[3]), int(linetemp[8])

    if e_1 - s_1 < 500:
        continue

    direction = linetemp[4]
 
    if id_1 != cid:
        if cid != '' and cid in nano_set:
            print_corrected_read(cid, seq_to_ranges, illu_to_ranges.setdefault(cid, []), length)
        cid = id_1
        length = len_1
        seq_to_ranges = {}

    if id_2 not in seq_to_ranges:
        seq_to_ranges[id_2] = (s_1, e_1, direction)
    else:
        str_1 = seq_to_ranges[id_2][0]
        str_2 = seq_to_ranges[id_2][1]
        d = seq_to_ranges[id_2][2]
        if direction == d and ( abs( str_1 - e_1 ) < 500 or abs(s_1 - str_2) < 500):
            seq_to_ranges[id_2] = (min(s_1, str_1), max(e_1, str_2), direction)

if cid in nano_set:    
    print_corrected_read(cid, seq_to_ranges, illu_to_ranges.setdefault(cid, []), length)
 
