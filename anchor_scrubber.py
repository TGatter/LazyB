#!/usr/bin/env python

import networkx as nx
import sys, os

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

paf_anchors = sys.argv[1] # anchors to nanopores
seq_file_base = sys.argv[2]
output_file = sys.argv[3]
tmp_dir = sys.argv[4]
output = open(output_file, "w")

if seq_file_base.endswith("fa") or seq_file_base.endswith("fasta"):
    seq_dict = SeqIO.index_db(seq_file_base+".idx", seq_file_base, "fasta")
else:
    seq_dict = SeqIO.index_db(seq_file_base+".idx", seq_file_base, "fastq")

def find_sequence_r(seqid, start, end):
   return seq_dict[seqid].seq[start:end + 1]

def find_sequence(seqid):
   return seq_dict[seqid].seq

def reverse_complement(seq):
   return seq.reverse_complement()

def write_sequence(f, seq):
    for i in range(0, len(seq), 60):
        f.write(str(seq[i:i + 60]) + "\n")


# Collect overlaps for all 
nano_to_ranges = {}
illu_to_nanos = {}

for line in open(paf_anchors, "r"):

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

    if id_2 not in nano_to_ranges:
        nano_to_ranges[id_2] = {}

    if id_1 not in nano_to_ranges[id_2]:
        nano_to_ranges[id_2][id_1] = (s_2, e_2)

    illu_to_nanos.setdefault(id_1, set()).add(id_2)

def print_corrected_read(cid, seq_to_ranges, nano_to_ranges, length):

    #print "Print", cid

    join = []
    for sid in seq_to_ranges:
        s, e, _ = seq_to_ranges[sid]
        join.append( (s,e) )
    for sid in nano_to_ranges:
        s, e = nano_to_ranges[sid]
        join.append( (s,e) )

    join.sort()

    #print cid, join

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
        write_sequence(output, find_sequence_r(cid, max(cs, 200), min(ce, length - 200)) )


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

for nano_id_chunk in chunker([ nid for nid in nano_to_ranges ], 30000) : # only do actual anchored nanos!

    other_nano_set = set()
    for nano_id in nano_id_chunk:
        for illu in nano_to_ranges[nano_id]:
            for nano in illu_to_nanos[illu]:
                other_nano_set.add(nano)

    # write out overlap group!
    f1 = open(tmp_dir+"/temp_reference.fa", "w")
    for nano_id in nano_id_chunk:
        f1.write(">"+nano_id+"\n")
        write_sequence(f1, find_sequence(nano_id))
    f1.close()

    f2 = open(tmp_dir+"/temp_others.fa", "w")
    for o_nano_id in other_nano_set:
        f2.write(">"+o_nano_id+"\n")
        write_sequence(f2, find_sequence(o_nano_id))
    f2.close()
    
    #print "A", nano_id_chunk
    #print "B", other_nano_set
    
    os.system('minimap2 -x ava-ont --dual=yes  '+tmp_dir+ '/temp_others.fa ' +tmp_dir+ '/temp_reference.fa > ' + tmp_dir+'/temp_pwa.paf 2> /dev/null')
 
    seq_to_ranges = {}
    length = 0
    cid = ""

    for line in open(tmp_dir+'/temp_pwa.paf', "r"):


        linetemp = line.rstrip().split("\t")
        if len(linetemp) == 1:
            continue;

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
            if cid != '':
                print_corrected_read(cid, seq_to_ranges, nano_to_ranges.setdefault(cid, []), length)
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
    
    print_corrected_read(cid, seq_to_ranges, nano_to_ranges.setdefault(cid, []), length)

