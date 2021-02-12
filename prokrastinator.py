#!/usr/bin/env python


#
#scr/k70san2/thomas/yeast_nanopore
#[thomas@vodka yeast_nanopore]$ ~/work/Prokrastinator/PROKRASTINATOR/prokrastinator_paf2.py paf assembly/seed_120/s288c-unitigs.l300 /scr/k70san2/thomas/yeast_nanopore/YeastStrainsStudy/fastqs/ont/s288c/data_Sanger_FAA86705_s288c_5_a_d_pass2D

# libraries
import networkx as nx
from networkx.utils import UnionFind
from networkx.algorithms.approximation import clique as nxc_a
from networkx.algorithms import clique as nxc_e

import collections
from collections import deque
import matplotlib.pyplot as plt # for drawing graphs
import itertools as it
import numpy as np
import sys
import operator
import math

from sortedcontainers import SortedSet

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

#from timeit import default_timer as timer

#import os
#import psutil

threshold_matches = 500
threshold_length = 500

blast_file = sys.argv[1]
illumina_file_base = sys.argv[2]
nano_file_base = sys.argv[3]
tmp_dir = sys.argv[4]

if len(sys.argv) > 5:
    threshold_matches = int(sys.argv[5])
if len(sys.argv) > 6:
    threshold_length = int(sys.argv[6])

min_matches = threshold_matches - threshold_matches * 0.2
min_length = threshold_length - threshold_length * 0.2

wiggle_room = 300
if len(sys.argv) > 7:
    wiggle_room = sys.argv[7]


output_base = tmp_dir+"/temp_1"
of_query = open(output_base + ".query.fa", "w")
of_paf = open(output_base + ".align.paf", "w")
of_target = open(output_base + ".target.fa", "w")

of_log = open(tmp_dir + "/prokrast_gen.log", "w")

# define node and edge datatypes

# illumina = source; nano = target
NodeMeta = collections.namedtuple('NodeMeta', 'nano_range illumina_range ratio direction score primary')
EdgeMeta = collections.namedtuple('EdgeMeta', 'overlap direction score primary')

Order = collections.namedtuple('Order', 'start end left_offset right_offset contained base_node score id_set direction primary')

ContainElement = collections.namedtuple('ContainElement', 'nano matches direction nano_length score primary')

# directed graph
DirInfo = collections.namedtuple('DirInfo', 'id nano_range_s illumina_range_s ratio_s overlap nano_range_e illumina_range_e ratio_e')

# flattened_regions
Region = collections.namedtuple('Region', 'a_start a_end id reverse is_nano s_start s_end')

# open up sequence map for illumina
illumina_dict = SeqIO.index_db(illumina_file_base+".idx", illumina_file_base, "fasta")

if nano_file_base.endswith("fa") or nano_file_base.endswith("fasta"):
    nano_dict = SeqIO.index_db(nano_file_base+".idx", nano_file_base, "fasta")
else:
    nano_dict = SeqIO.index_db(nano_file_base+".idx", nano_file_base, "fastq")

def find_sequence_illumina(seqid, start, end):
   return illumina_dict[seqid].seq[start:end + 1]

def find_sequence_nano(seqid, start, end):
   return nano_dict[seqid].seq[start:end + 1]

def reverse_complement(seq):
   return seq.reverse_complement()

def write_sequence(f, seq):
    for i in range(0, len(seq), 60):
        f.write(str(seq[i:i + 60]) + "\n")

# create graph

G = nx.Graph()

prevHitId = ""  # just ot prevent the error at the first line
chunkNodes = {}

node_index = 0  # we give each node a global index we can write out!

id_to_overlap = {}

#process = psutil.Process(os.getpid())
#print "BASE", process.memory_info().rss  # in bytes 

############## Read Blast in and build base graph

#start = timer()

for line in open(blast_file, "r"):

    # read line into a list

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

    linetemp = line.rstrip().split("\t")
    if len(linetemp) == 1:
        continue;

    # get query and subject IDs
    illuminaId, nanoId = linetemp[0], linetemp[5]  # remember about 0-based counting
    # get coordinates in both sequences + score
    illuminaRange = ( int(linetemp[2]) , int(linetemp[3]) - 1 )  # remember about 0-based to 1-based counting
    nanoRange = ( int(linetemp[7]) , int(linetemp[8]) - 1 )  # order in coordinates defines hit direction

    if linetemp[4] == '-':
        direct = -1  # we use an actual int to make life a little easier as we can use the normal *
    else:
        direct = +1

    matches = int(float(linetemp[9])) # number of actual matches over the length

    nanoLength = int(linetemp[6])

    rratio = (illuminaRange[1] - illuminaRange[0] + 1) / float(nanoRange[1] - nanoRange[0] + 1)  # illumina tends to be longer

    # compatible?
    if matches < min_matches:
        continue
    if illuminaRange[1] - illuminaRange[0] + 1 <  min_length:
        continue

    primary = illuminaRange[1] - illuminaRange[0] + 1 >=  threshold_length and matches >= threshold_matches

    # add node to the big graph.
    # Remember: the table needs to be sorted by scaffold!
    # if this node already exists, then add attributes to the list
    if not G.has_node(nanoId):  # if this node does not exist, then create it
        G.add_node(nanoId, print_index=node_index)
        G.nodes[nanoId]['nano_length'] = nanoLength # this value should be conistent in input between lines!
        node_index = node_index + 1  # each node gets a unique (small!) id as well
        G.nodes[nanoId]['matches'] = {} # initialize empty matches

    # append the new inf
    if illuminaId in G.nodes[nanoId]['matches']:
        continue  # TODO: for now kill double hits

    G.nodes[nanoId]['matches'][illuminaId] = NodeMeta(nano_range=nanoRange, illumina_range=illuminaRange, ratio=rratio, direction=direct, score=matches, primary=primary)

    # are we still in the same scaffold? If not, let's start a new one
    if illuminaId != prevHitId:
        chunkNodes = []
        prevHitId = illuminaId

    # nano_id is a unique identifier, and we can retrieve nodes from the graph in O(1) with it

    # here are reads associated with the same contig
    # and now iterate only through these reads
    for prevId in chunkNodes:

        newQ = G.nodes[nanoId]['matches'][illuminaId].illumina_range
        prevQ = G.nodes[prevId]['matches'][illuminaId].illumina_range

        # the interval of the query both share
        overlapQ = (max(newQ[0], prevQ[0]), min(newQ[1], prevQ[1]))

        is_end_a =  (G.nodes[nanoId]['nano_length'] - G.nodes[nanoId]['matches'][illuminaId].nano_range[1] ) < 600
        is_start_a =  G.nodes[nanoId]['matches'][illuminaId].nano_range[0] < 600

        if G.nodes[nanoId]['matches'][illuminaId].direction < 0:
            tmp = is_end_a
            is_end_a = is_start_a
            is_start_a = tmp

        is_end_b =  (G.nodes[prevId]['nano_length'] - G.nodes[prevId]['matches'][illuminaId].nano_range[1] ) < 600
        is_start_b =  G.nodes[prevId]['matches'][illuminaId].nano_range[0] < 600

        if G.nodes[prevId]['matches'][illuminaId].direction < 0:
            tmp = is_end_b
            is_end_b = is_start_b
            is_start_b = tmp

        bad = False
        #bad = ( (overlapQ[0] - newQ[0] > 500 or overlapQ[0] - prevQ[0] > 500 ) and not is_start_a and not is_start_b ) or ( (prevQ[1] - overlapQ[1] > 500 or newQ[1] - overlapQ[1] > 500) and not is_end_a and not is_end_b )
        #if bad and overlapQ[0] <= overlapQ[1] and overlapQ[1] - overlapQ[0] > 100:
        #    print nanoId, prevId, illuminaId, bad
        #    print nanoId, newQ, G.nodes[nanoId]['matches'][illuminaId].nano_range, G.nodes[nanoId]['nano_length'], is_start_a, is_end_a, G.nodes[nanoId]['matches'][illuminaId].direction
        #    print prevId, prevQ, G.nodes[prevId]['matches'][illuminaId].nano_range, G.nodes[prevId]['nano_length'], is_start_b, is_end_b, G.nodes[prevId]['matches'][illuminaId].direction

        # and if the alignment coordinates within the query contig intersect:
        if overlapQ[0] <= overlapQ[1] and overlapQ[1] - overlapQ[0] > 100 and not bad:  # TODO: should we care for the size of the interval here ?
            # then get the relative direction
            edge_direct = G.nodes[nanoId]['matches'][illuminaId].direction * G.nodes[prevId]['matches'][illuminaId].direction
            is_primary = G.nodes[nanoId]['matches'][illuminaId].primary * G.nodes[prevId]['matches'][illuminaId].primary

            # and lengths of alignment for each read
            thisLen = newQ[1] - newQ[0] + 1
            prevLen = prevQ[1] - prevQ[0] + 1
            common_length = overlapQ[1] - overlapQ[0] + 1

            # NP overlap score # TODO Attention, this deliberately rounds to int!
            thisScore = G.nodes[nanoId]['matches'][illuminaId].score * common_length / thisLen
            prevScore = G.nodes[prevId]['matches'][illuminaId].score * common_length / prevLen

            overlap_score_min = min(thisScore, prevScore)
            overlap_score_sum = thisScore + prevScore

            # add edge if necessary
            if not G.has_edge(prevId, nanoId):
                G.add_edge(prevId, nanoId)
                G.edges[prevId, nanoId]['matches'] = {}

            G.edges[prevId, nanoId]['matches'][illuminaId] = EdgeMeta(overlap=overlapQ, direction=edge_direct, score=overlap_score_sum, primary=is_primary)

        # we add this after the loop to avoid the identity check
    chunkNodes.append(nanoId)  # this is a unique identifier, and we can retrieve nodes from the graph in O(1) with it

of_log.write("# Primary Nodes " + str(len(G.nodes())) + ".\n")
of_log.write("# Primary Edges " + str(len(G.edges())) + ".\n")

# clean up
chunkNodes = {}
prevHitId = ""

#print "Reading PAF", process.memory_info().rss  # in bytes 

#end = timer()
#print "Reading PAF", (end - start)

#for n in G.nodes():
#    print "NodeMatch", n,
#    for match_id in G.nodes[n]['matches']:
#        print match_id,
#    print

############## Check simple compatibility within edges

def matches_compatible(G, edge, il_id_1, il_id_2):

    ##### check distance

    d1 = G.edges[edge]['matches'][il_id_1].direction
    d2 = G.edges[edge]['matches'][il_id_2].direction

    # shared data

    overlap_1 = G.edges[edge]['matches'][il_id_1].overlap
    overlap_2 = G.edges[edge]['matches'][il_id_2].overlap

    # compute on nano 1

    ratio_1_1 = G.nodes[edge[0]]['matches'][il_id_1].ratio
    ratio_1_2 = G.nodes[edge[0]]['matches'][il_id_2].ratio

    nc_1_1_l = (overlap_1[0] -  G.nodes[edge[0]]['matches'][il_id_1].illumina_range[0]) / ratio_1_1
    nc_1_1_r = (G.nodes[edge[0]]['matches'][il_id_1].illumina_range[1] - overlap_1[1]) / ratio_1_1

    if G.nodes[edge[0]]['matches'][il_id_1].direction < 0: # flip
        t = nc_1_1_l
        nc_1_1_l = nc_1_1_r
        nc_1_1_r = t

    nc_1_2_l = (overlap_2[0] -  G.nodes[edge[0]]['matches'][il_id_2].illumina_range[0]) / ratio_1_2
    nc_1_2_r = (G.nodes[edge[0]]['matches'][il_id_2].illumina_range[1] - overlap_2[1]) / ratio_1_2

    if G.nodes[edge[0]]['matches'][il_id_2].direction < 0: # flip
        t = nc_1_2_l
        nc_1_2_l = nc_1_2_r
        nc_1_2_r = t

    nano_1_1 = G.nodes[edge[0]]['matches'][il_id_1].nano_range
    nano_1_2 = G.nodes[edge[0]]['matches'][il_id_2].nano_range
    corrected_nano_1_1 = ( G.nodes[edge[0]]['matches'][il_id_1].nano_range[0] + nc_1_1_l, G.nodes[edge[0]]['matches'][il_id_1].nano_range[1] - nc_1_1_r)
    corrected_nano_1_2 = ( G.nodes[edge[0]]['matches'][il_id_2].nano_range[0] + nc_1_2_l, G.nodes[edge[0]]['matches'][il_id_2].nano_range[1] - nc_1_2_r)

    if corrected_nano_1_1[0] <= corrected_nano_1_2[1] and corrected_nano_1_2[0] <= corrected_nano_1_1[1] : # regions are overlapping!
        orientation_1  = 0
        diff_1 = 0
        # check kind of overlap
        if corrected_nano_1_1[0] < corrected_nano_1_2[0] and corrected_nano_1_1[1] < corrected_nano_1_2[1]: # 1 falls naturally left of 2   
            orientation_1  = 2
            diff_1 = corrected_nano_1_1[1] - corrected_nano_1_2[0] + 1
        elif corrected_nano_1_1[0] > corrected_nano_1_2[0] and corrected_nano_1_1[1] > corrected_nano_1_2[1]: # 2 falls naturally left of 1
            orientation_1  = -2
            diff_1 = corrected_nano_1_2[1] - corrected_nano_1_1[0] + 1
    elif corrected_nano_1_1[0] < corrected_nano_1_2[0]: # ID 1 naturally left of 2
        orientation_1  = 1
        diff_1 = corrected_nano_1_2[0] - corrected_nano_1_1[1] + 1
    else: # ID 2 naturally left of 2
        orientation_1  = -1
        diff_1 = corrected_nano_1_1[0] - corrected_nano_1_2[1] + 1

    if nano_1_1[0] <= nano_1_2[1] and nano_1_2[0] <= nano_1_1[1] : # regions are overlapping!
        if nano_1_1[0] < nano_1_2[0] and nano_1_1[1] < nano_1_2[1]: # 1 falls naturally left of 2   
            uco = 2
        elif nano_1_1[0] > nano_1_2[0] and nano_1_1[1] > nano_1_2[1]: # 2 falls naturally left of 1
            uco = -2
        else:
            uco = 0

        if orientation_1 < 0 and not uco < 0:
            return False
        elif orientation_1 > 0 and not uco > 0:
            return False

    # same on nano 2

    ratio_2_1 = G.nodes[edge[1]]['matches'][il_id_1].ratio
    ratio_2_2 = G.nodes[edge[1]]['matches'][il_id_2].ratio

    nc_2_1_l = (overlap_1[0] -  G.nodes[edge[1]]['matches'][il_id_1].illumina_range[0]) / ratio_2_1
    nc_2_1_r = (G.nodes[edge[1]]['matches'][il_id_1].illumina_range[1] - overlap_1[1]) / ratio_2_1

    if G.nodes[edge[1]]['matches'][il_id_1].direction < 0: # flip
        t = nc_2_1_l
        nc_2_1_l = nc_2_1_r
        nc_2_1_r = t

    nc_2_2_l = (overlap_2[0] -  G.nodes[edge[1]]['matches'][il_id_2].illumina_range[0]) / ratio_2_2
    nc_2_2_r = (G.nodes[edge[1]]['matches'][il_id_2].illumina_range[1] - overlap_2[1]) / ratio_2_2

    if G.nodes[edge[1]]['matches'][il_id_2].direction < 0: # flip
        t = nc_2_2_l
        nc_2_2_l = nc_2_2_r
        nc_2_2_r = t

    nano_2_1 = G.nodes[edge[1]]['matches'][il_id_1].nano_range
    nano_2_2 = G.nodes[edge[1]]['matches'][il_id_2].nano_range
    corrected_nano_2_1 = ( G.nodes[edge[1]]['matches'][il_id_1].nano_range[0] + nc_2_1_l, G.nodes[edge[1]]['matches'][il_id_1].nano_range[1] - nc_2_1_r )
    corrected_nano_2_2 = ( G.nodes[edge[1]]['matches'][il_id_2].nano_range[0] + nc_2_2_l, G.nodes[edge[1]]['matches'][il_id_2].nano_range[1] - nc_2_2_r )

    if corrected_nano_2_1[0] <= corrected_nano_2_2[1] and corrected_nano_2_2[0] <= corrected_nano_2_1[1] : # regions are overlapping!
        orientation_2  = 0
        diff_2 = 0
        if corrected_nano_2_1[0] < corrected_nano_2_2[0] and corrected_nano_2_1[1] < corrected_nano_2_2[1]: # 1 falls naturally left of 2   
            orientation_2  = 2
            diff_2 = corrected_nano_2_1[1] - corrected_nano_2_2[0] + 1
        elif corrected_nano_2_1[0] > corrected_nano_2_2[0] and corrected_nano_2_1[1] > corrected_nano_2_2[1]: # 2 falls naturally left of 1
            orientation_2  = -2
            diff_2 = corrected_nano_2_2[1] - corrected_nano_2_1[0] + 1
    elif corrected_nano_2_1[0] < corrected_nano_2_2[0]: # ID 1 naturally left of 2
        orientation_2  = 1
        diff_2 = corrected_nano_2_2[0] - corrected_nano_2_1[1] + 1
    else: # ID 2 naturally left of 2
        orientation_2  = -1
        diff_2 = corrected_nano_2_1[0] - corrected_nano_2_2[1] + 1

    if nano_2_1[0] <= nano_2_2[1] and nano_2_2[0] <= nano_2_1[1] : # regions are overlapping!
        if nano_2_1[0] < nano_2_2[0] and nano_2_1[1] < nano_2_2[1]: # 1 falls naturally left of 2   
            uco = 2
        elif nano_2_1[0] > nano_2_2[0] and nano_2_1[1] > nano_2_2[1]: # 2 falls naturally left of 1
            uco = -2
        else:
            uco = 0

        if orientation_2 < 0 and not uco < 0:
            return False
        elif orientation_2 > 0 and not uco > 0:
            return False

    # now see if they are compatible

    if d1 < 0: # we arbitrarily switch the second
        orientation_2 = orientation_2 * -1

    if orientation_1 == orientation_2 and orientation_1 != 0: # easy case, just check for the
        mm_diff = max(diff_1, diff_2) - min(diff_1, diff_2)
        matching = ( mm_diff <= wiggle_room ) or ( mm_diff * 100 / max(diff_1, diff_2) <= 15 )  # TODO > 15 % seems to be a decent cut, could be < 20% as well
    elif ( orientation_1 < 0 and orientation_2 < 0 ) or ( orientation_1 > 0 and orientation_2 > 0 ):
        matching = diff_1 + diff_2 <= wiggle_room
    else:
        matching = False

    return matching

def max_pairwise_paths(G, edge, ids, direction):

    if len(ids) == 0:
        return []

    #print "max_pairwise_paths", edge

    s_list = []
    e_list = []
    for il_id in ids:       # we gather up all matches, before searching for stretches
        s = G.nodes[edge[0]]['matches'][il_id]
        e = G.nodes[edge[1]]['matches'][il_id]
        s_list.append( (s.nano_range, il_id) )
        e_list.append( (e.nano_range, il_id) )

    # sort them by position on nanopores
    s_list.sort()
    e_list.sort()
    if direction < 0: # arbitrarily switch the second order
        e_list.reverse()

    # we do some dynamic programming here (size must be the same)

    length = len(s_list)
    population = [ ([], G.edges[edge]['matches'][x].score) for _ , x in s_list ]
    for k in range(length - 1):
        for l in range(k + 1, length):

            compatible = matches_compatible(G, edge, s_list[k][1], s_list[l][1])
            score = population[k][1] +  G.edges[edge]['matches'][s_list[l][1]].score

            if compatible and  score > population[l][1]:
                population[l] = ( population[k][0] + [k], score )

    max_path_val = 0
    max_path_i = 0
    for k in range(len(population)):
        population[k][0].append(k)   # add in last index while finding max score
        if population[k][1] > max_path_val:
            max_path_val = population[k][1]
            max_path_i = k

    viable_list = []
    viable_list.append( ([ s_list[i][1] for i in population[max_path_i][0] ] , max_path_val, any(G.edges[edge]['matches'][s_list[i][1]].primary for i in population[max_path_i][0]) or len(population[max_path_i][0]) > 2 ) )

    for k in range(len(population)):
        if population[k][1] > max_path_val * 0.75:  # within 75% of score to optimal
            disjoint = True
            for vl in viable_list:
                disjoint = disjoint and not any(s_list[i][1] in vl[0] for i in population[k][0])
            if disjoint:
                n_ids = [ s_list[i][1] for i in population[k][0] ]
                viable_list.append( (n_ids, population[k][1], any(G.edges[edge]['matches'][s_list[i][1]].primary for i in population[k][0]) ) )

    if len(viable_list) == 1 and viable_list[0][2] == True: # if other, definitively already a shadow

        s_ids = []
        for match_id in G.nodes[edge[0]]['matches']:
            s_ids.append( ( G.nodes[edge[0]]['matches'][match_id].nano_range, match_id ) )
        e_ids = []
        for match_id in G.nodes[edge[1]]['matches']:
            e_ids.append( ( G.nodes[edge[1]]['matches'][match_id].nano_range, match_id ) )

        s_ids.sort()
        e_ids.sort()

        #print "raw ", s_ids, e_ids

        s_ids_ = s_ids # filter containments
        #for r, i in s_ids:
        #    if len(s_ids_) == 0 or s_ids_[-1][0][1] < r[1]:
        #        s_ids_.append( (r, i) )

        e_ids_ = e_ids # filter containments
        #for r, i in e_ids:
        #    if len(e_ids_) == 0 or e_ids_[-1][0][1] < r[1]:
        #        e_ids_.append( (r, i) )

        if direction < 0: # arbitrarily switch the second order, consistent with above
            e_ids_.reverse()

        ids = viable_list[0][0]
        if (s_ids_[0][1] != ids[0] and e_ids_[0][1] != ids[0]) or (s_ids_[-1][1] != ids[-1] and e_ids_[-1][1] != ids[-1]): # diverging starts / ends
            #print "BAD A", ids, s_ids_, e_ids_
            return [ (viable_list[0][0], viable_list[0][1], False) ]
        else:
            #print "Ch", ids, s_ids_, e_ids_
            i = 0
            j = 0
            shadow = False
            for id_ in ids:
                inter_a = False
                while id_ != s_ids_[i][1]:
                    #print "A", i, id_, s_ids_[i][1]
                    i = i + 1
                    inter_a = True
                i = i + 1

                inter_b = False
                while id_ != e_ids_[j][1]:
                    #print "B", j, id_, e_ids_[j][1]
                    j = j +1
                    inter_b = True
                j = j + 1
                if inter_a and inter_b:
                    shadow = True
                    break
            if shadow:
                #print "BAD B"
                return [ (viable_list[0][0], viable_list[0][1], False) ]

    return viable_list


############## Compute Overlaps

def compute_overhangs(G, il_id, edge, node):

    nano_len = G.nodes[node]['nano_length']
    nano_range = G.nodes[node]['matches'][il_id].nano_range

    nano_correction_l = ( G.edges[edge]['matches'][il_id].overlap[0] - G.nodes[node]['matches'][il_id].illumina_range[0] ) / G.nodes[node]['matches'][il_id].ratio
    nano_correction_r = ( G.nodes[node]['matches'][il_id].illumina_range[1] - G.edges[edge]['matches'][il_id].overlap[1] ) / G.nodes[node]['matches'][il_id].ratio

    if G.nodes[node]['matches'][il_id].direction < 0:
        t = nano_correction_l
        nano_correction_l = nano_correction_r
        nano_correction_r = t


    left_overhang = nano_range[0] + nano_correction_l
    right_overhang = nano_len - nano_range[1] + nano_correction_r

    #if G.nodes[node]['matches'][il_id].direction < 0:
    #    t = left_overhang
    #    left_overhang = right_overhang
    #    right_overhang = t
    
    return left_overhang, right_overhang


def compute_overlap_from_id_set(G, id_set, edge, direction, score, primary):

    id_l = id_set[0]
    id_r = id_set[-1]

    lo_l_1, ro_l_1  = compute_overhangs(G, id_l, edge, edge[0])  
    lo_r_1, ro_r_1 = compute_overhangs(G, id_r, edge, edge[0])

    lo_l_2, ro_l_2 = compute_overhangs(G, id_l, edge, edge[1])
    lo_r_2, ro_r_2 = compute_overhangs(G, id_r, edge, edge[1])

    left_overhang_1 = lo_l_1
    right_overhang_1 = ro_r_1

    if direction < 0:
        left_overhang_2 = ro_l_2
        right_overhang_2 = lo_r_2
    else:
        left_overhang_2 = lo_l_2
        right_overhang_2 = ro_r_2

    # check direction implied, as we corrected the matched regions, they cover the same stretch, hence overhangs can be directly compared
    if left_overhang_1 <= left_overhang_2 and right_overhang_1 <= right_overhang_2:   # for overlaps always from contained to container, setting edge[0] to 1+ for l_id ! 
        return Order(start=edge[0] , end=edge[1] , left_offset=left_overhang_2-left_overhang_1, right_offset=right_overhang_2-right_overhang_1, contained = True, base_node=edge[0], score = score, id_set = id_set, direction = direction, primary = primary)
    elif left_overhang_1 >= left_overhang_2 and right_overhang_1 >= right_overhang_2:
        return Order(start=edge[1] , end=edge[0] , left_offset=left_overhang_1-left_overhang_2, right_offset=right_overhang_1-right_overhang_2, contained = True, base_node=edge[0], score = score, id_set = id_set, direction = direction, primary = primary)
    elif left_overhang_1 > left_overhang_2 and right_overhang_1 < right_overhang_2:  # for overhang we go from perceived left to right, setting edge[0] as +1 for l_id!
        return Order(start=edge[0] , end=edge[1] , left_offset=left_overhang_1-left_overhang_2, right_offset=right_overhang_2-right_overhang_1, contained = False, base_node=edge[0], score = score, id_set = id_set, direction = direction, primary = primary)
    elif left_overhang_1 < left_overhang_2 and right_overhang_1 > right_overhang_2:
        return Order(start=edge[1] , end=edge[0] , left_offset=left_overhang_2-left_overhang_1, right_offset=right_overhang_1-right_overhang_2, contained = False, base_node=edge[0], score = score, id_set = id_set, direction = direction, primary = primary)

#start = timer()

for edge in G.edges():

    # we aim to find the longest stretch of compatible matches

    # find direction-compatible sets
    minus_ids = []
    plus_ids = []
    for il_id in G.edges[edge]['matches']:
        if G.edges[edge]['matches'][il_id].direction < 0:
            minus_ids.append(il_id)
        else:
            plus_ids.append(il_id)

    minus_paths = max_pairwise_paths(G, edge, minus_ids, -1) # find the longest compatible path for this one
    plus_paths = max_pairwise_paths(G, edge, plus_ids, 1) # find the longest compatible path for this one
    # [ ID_List, bit_score, primary]

    # first translate to     
    has_primary = any( primary for _ , _, primary in minus_paths ) or any( primary for _ , _, primary in plus_paths )
    if has_primary:  # filter out only primaries
        minus_paths[:] = [ e for e in minus_paths if e[2]]
        plus_paths[:] = [ e for e in plus_paths if e[2]]

    has_multi = any( len(ids) > 1 for ids , _, _ in minus_paths ) or any( len(ids) > 1 for ids , _, _ in plus_paths )
    if has_multi:  # filter out only multis
        minus_paths[:] = [ e for e in minus_paths if len(e[0]) > 1]
        plus_paths[:] = [ e for e in plus_paths if len(e[0]) > 1]    

    # everything that can be safely deleted is, so now test if we need to shadow ban this edge!
    if len(minus_paths) + len(plus_paths) > 1: # we have at least 2 (possibly primary) sets, definitely something wrong
        G.edges[edge]['shadow'] = True
        # for now we keep all
        #if len(minus_paths) > 0 and has_multi: # no saming on same direction if we have multis!
        #    # we keep only the max one and block it as secondary
        #    max_id = max(minus_paths, key=lambda item: item[1])
        #    max_id[2] = False
        #    minus_paths = [max_id] # reset to best
        #if len(plus_paths) > 0 and has_multi: # no saming on same direction if we have multis!
        #    # we keep only the max one and block it as secondary
        #    max_id = max(plus_paths, key=lambda item: item[1])
        #    max_id[2] = False
        #    plus_paths = [max_id] # reset to best
    else: 
        path = (minus_paths + plus_paths)[0]
        G.edges[edge]['shadow'] = not path[2]
    
    G.edges[edge]['orders'] = []
    for path in minus_paths:
        G.edges[edge]['orders'].append( compute_overlap_from_id_set(G, path[0], edge, -1, path[1], path[2]) )

    for path in plus_paths:
        G.edges[edge]['orders'].append( compute_overlap_from_id_set(G, path[0], edge, 1,  path[1], path[2]) )

#print "Chains and Overlaps", process.memory_info().rss  # in bytes 
#end = timer()
#print "Chains and Overlaps", (end - start)

############## We Contract Contains, correct, and unify orders

def sanity_check_matching_edge(G, subnode, node, target, contain_order):

    contain_edge = (node, subnode)
    checkon = (node, target)
    checkfor = (subnode, target)

    for i, checkon_order in enumerate(G.edges[checkon]['orders']):
        for j, checkfor_order in enumerate(G.edges[checkfor]['orders']):

            sane = True
            sane = sane and contain_order.direction * checkon_order.direction == checkfor_order.direction  # directions match up

            c1 = checkfor_order.contained
            c2 = checkon_order.contained

            if c1 and c2: # both contained
                sane = sane and ( (checkfor_order.end == target or checkfor_order.start == target) and checkon_order.start == target)
            elif c1 and not c2:

                if checkfor_order.end != target: # otherwise always OK!

                    if ( checkon_order.direction < 0 and ( (node == checkon_order.base_node and checkon_order.end == target) or (node != checkon_order.base_node and checkon_order.end == target)) ) \
                          or (checkon_order.direction > 0 and checkon_order.end == target):

                        if contain_order.direction < 0:
                            l1 = False     # contain_order
                            l2 = True      # checkfor_order
                            l3 = False     # checkon_order
                        else:
                            l1 = False
                            l2 = False
                            l3 = False
                    else:
                        if contain_order.direction < 0:
                            l1 = True
                            l2 = False
                            l3 = True
                        else:
                            l1 = True
                            l2 = True
                            l3 = True

                    if contain_order.direction < 0 and contain_order.base_node != contain_order.end:
                        l1 = not l1
                    if checkfor_order.direction < 0 and checkfor_order.base_node != checkfor_order.end:
                        l2 = not l2

                    d1 = contain_order.left_offset if l1 else contain_order.right_offset
                    d2 = checkfor_order.left_offset if l2 else checkfor_order.right_offset
                    d3 = checkon_order.left_offset if l3 else checkon_order.right_offset

                    sane = sane and ( d1+d2+d3 < wiggle_room)

            elif not c1 and c2:
                sane = sane and (checkon_order.start == target)
            else:  # both are not contained; overlap, check

                d1 = checkfor_order.start == target
                d2 = checkon_order.start == target

                if checkfor_order.direction < 0 and checkfor_order.base_node == target: # flip if needed
                    d1 = not d1

                if checkon_order.direction < 0 and checkon_order.base_node == target:
                    d2 = not d2

                if contain_order.direction < 0: # (double) flip sub-edge
                    d1 = not d1

                sane = sane and d1 == d2
            if sane:
                return (True, i, j)

    return (False, None, None)

#start = timer()

# find edges to contract
contraction_edges = {}
for edge in G.edges():
    for order in G.edges[edge]['orders']:
        if order.contained and order.primary : # we have a contained edge, do the contraction with sanity checks! (we contract sane shadows)

            subnode = order.start # the contained node
            node = order.end # the container node

            sane = True
            sanity_set = []
            for checkfor in G.edges(subnode):

                target = checkfor[1]

                if target == node or G.edges[checkfor]['shadow']:
                    continue

                sane = sane and G.has_edge(node, target)
                if not sane:
                    break

                s, i, j = sanity_check_matching_edge(G, subnode, node, target, order)
                sane = sane and s

                if not sane:
                    break

                sanity_set.append( (checkfor, j, (node, target), i) ) # for, order, on, order 

            if sane:
                contraction_edges[edge] = (sanity_set, order)
                break



#for cc in nx.connected_component_subgraphs(G):
#  print("RAW ------------------------------------------------")
#  for edge in cc.edges(): 
#     print edge, edge in contraction_edges or (edge[1], edge[0]) in contraction_edges, cc.edges[edge]

# mark nodes to contract into
for node in G.nodes():
    G.nodes[node]['contract_to'] = node

for contain in contraction_edges:
    sanity_set, contain_order = contraction_edges[contain]

    contract_to = G.nodes[contain_order.end]['contract_to']  #  get which node we contract this one to, set recursively!

    if G.nodes[contain_order.start]['contract_to'] == contain_order.start:  # this node was untouched so far (could be contracted into!)
        G.nodes[contain_order.start]['contract_to'] = contract_to  # set contraction for next time

delete_nodes = [] 
for contain in contraction_edges:
    sanity_set, contain_order = contraction_edges[contain]

    delete_nodes.append(contain_order.start) # mark for later delete

    contract_to = G.nodes[contain_order.start]['contract_to']  # move the root token
    G.nodes[contract_to]['contract_root'] = True

    if 'contract_root' in G.nodes[contain_order.start]:
       del G.nodes[contain_order.start]['contract_root']
    

# do the contraction
for node in G.nodes():  # set empty lists for all elements
    G.nodes[node]['contain_elements'] = []

for contain in contraction_edges:
    sanity_set, contain_order = contraction_edges[contain]

    if 'contract_root' not in G.nodes[contain_order.end]:
        continue 

    matches = {}
    for id_i in contain_order.id_set:
        matches[id_i] = G.nodes[contain_order.start]['matches'][id_i] # we don't actually need the edge overlaps, as we have anchors for that

    G.nodes[contain_order.end]['contain_elements'].append(ContainElement(nano = contain_order.start, matches = matches, direction = contain_order.direction, nano_length = G.nodes[contain_order.start]['nano_length'], score = contain_order.score, primary = contain_order.primary))

    # contain edge is always unique

for node in delete_nodes:
    if G.has_node(node):
        G.remove_node(node)

of_log.write("# Nodes after Collapsed Contains " + str(len(G.nodes())) + ".\n")
of_log.write("# Edges after Collapsed Contains " + str(len(G.edges())) + ".\n")

#for cc in nx.connected_component_subgraphs(G):
#  print("Contracted  ------------------------------------------------ ")
#  for edge in cc.edges(): 
#     print edge, cc.edges[edge], 
#  print("CONTAINS  -------- ")
#  for n in cc.nodes():
#     print "C", n, cc.nodes[n]['contain_elements']

############## Remove Edges with Contains

delete_edges = []
for edge in G.edges():

    G.edges[edge]['orders'][:] = [o for o in G.edges[edge]['orders'] if not o.contained]

    empty = len(G.edges[edge]['orders']) == 0
    if empty:
        delete_edges.append(edge)

for edge in delete_edges:
    G.remove_edge(edge[0], edge[1])

of_log.write("# Nodes after Bad Contains " + str(len(G.nodes())) + ".\n")
of_log.write("# Edges after Bad Contains " + str(len(G.edges())) + ".\n")

#end = timer()
#print "Contracting", (end - start)

#for cc in nx.connected_component_subgraphs(G):
#  print("Contain_Removed ------------------------------------------------ ")
#  for edge in cc.edges(): 
#     print edge, cc.edges[edge]

############## Check for Obvious Minus Cycles on non-showed regions!

#print "Contracting", process.memory_info().rss
#start = timer()

# compute the bitweight first!
shadow_count = 0
for edge in G.edges():
    G.edges[edge]['weight'] = 0
    if not G.edges[edge]['shadow']: # if not shadowed, we can max have one order!
        G.edges[edge]['weight'] = G.edges[edge]['orders'][0].score
        G.edges[edge]['consensus_direction'] = G.edges[edge]['orders'][0].direction
    else:
        shadow_count = shadow_count + 1
        if len(set( [ o.direction for o in G.edges[edge]['orders'] ] )) == 1:
            G.edges[edge]['consensus_direction'] = G.edges[edge]['orders'][0].direction

of_log.write("Shadow Count " + str(shadow_count) + ".\n")

def custom_kruskal_mst(G, ignore_nan=False):

    subtrees = UnionFind()
    edges = G.edges(data=True)

    def filter_nan_edges(edges=edges):
        sign = -1
        for u, v, d in edges:
            wt = d.get('weight', 1) * sign
            if np.isnan(wt) or 'consensus_direction' not in d:  # skip shadows!
                continue
            yield wt, u, v, d

    edges = sorted(filter_nan_edges(), key=operator.itemgetter(0))
    for wt, u, v, d in edges:
        if subtrees[u] != subtrees[v]:
            yield (u, v)
            subtrees.union(u, v)

def maximum_spanning_tree(G):

    edges = custom_kruskal_mst(G)
    edges = list(edges)
    T = G.__class__()  # Same graph class as G
    T.graph.update(G.graph)
    T.add_nodes_from(G.nodes.items())
    T.add_edges_from(edges)
    return T

# now Kruskal on graph without shadow
mst = maximum_spanning_tree(G)

remove_edges = []
inconsistent_minus_cycles = 0 
for edge in G.edges():
    if 'consensus_direction' in G.edges[edge] and not edge in mst.edges():
        path = nx.shortest_path(mst, edge[0], edge[1])
        direction = G.edges[edge[0], edge[1]]['consensus_direction']
        base_weight = G.edges[edge[0], edge[1]]['weight']
        weights = []
        for v, u in zip(path[:-1], path[1:]):
            direction = direction * G.edges[v,u]['consensus_direction']
            weights.append(G.edges[v,u]['weight']);
        if direction == -1:
            index_min = np.argmin(weights)
            weight_min = weights[index_min]
            index_max = np.argmax(weights)
            weight_max = weights[index_max]
            if weight_min < base_weight or ( base_weight * 1.1 >= weight_min and weight_min < weight_max * 0.8 ): #TODO double delete within 10% window if more than 20% less then max element
                remove_edges.append((path[index_min], path[index_min+1]))
                #print "Inconsistent Edge", (path[index_min], path[index_min+1])
            remove_edges.append(edge) # always remove the bad cycle edge
            inconsistent_minus_cycles += 1
        else:
            G.edges[edge[0], edge[1]]['cycle'] = 0

for edge in remove_edges:
    if G.has_edge(edge[0], edge[1]):
        G.remove_edge(edge[0], edge[1])

of_log.write("# Nodes Decycled " + str(len(G.nodes())) + ".\n")
of_log.write("# Edges Decycled " + str(len(G.edges())) + ".\n")

#end = timer()
#print "Cycle-Cancel", (end - start)

############## Create Directed Graph with possibly opposing overlap clusters

####### Direction Recursion

def direction_graph(cc, first_node, dg): # we do this as a normal DFS

    stack = deque()
    stack.append( (first_node, 1) )

    while len(stack) > 0:

        node, direction_mod = stack.pop()

        if not dg.has_node(node):
            dg.add_node(node) # add new node to digraph, and mark it as visited thus
            dg.nodes[node]['nano_length'] = cc.nodes[node]['nano_length'] 
            dg.nodes[node]['matches'] = {}
            dg.nodes[node]['contain_elements'] = G.nodes[node]['contain_elements']

        if 'direction' not in dg.nodes[node]:
            dg.nodes[node]['direction'] = direction_mod

        for edge in cc.edges(node):

            other = edge[1]   #other is always second

            other_exists = dg.has_node(other) and 'direction' in dg.nodes[other] # check before we add in any nodes

            # add in both nodes if needed
            if not dg.has_node(other):
                dg.add_node(other) # add new node to digraph, and mark it as visited thus 
                dg.nodes[other]['nano_length'] = cc.nodes[other]['nano_length']
                dg.nodes[other]['matches'] = {}
                dg.nodes[other]['contain_elements'] = G.nodes[other]['contain_elements']

            if dg.has_edge(edge[0], edge[1]) or dg.has_edge(edge[1], edge[0]): # this means we already added this cc edge in reverse
                continue

            for order in cc.edges[edge]['orders']:
                # check and set directions
                flip = False
                if order.direction < 0 and order.base_node == other:
                    flip = not flip
                if direction_mod < 0:
                    flip = not flip

                if flip == True:
                    s = order.end
                    e = order.start 
                else:
                    s = order.start
                    e = order.end

                add_meta = False
                if not dg.has_edge(s,e):
                    dg.add_edge(s , e)
                    dg.edges[s,e]['matches'] = {}
                    dg.edges[s,e]['orders'] = []
                    dg.edges[s,e]['shadow'] = G.edges[edge]['shadow']
                    if 'cycle' in  G.edges[edge]:
                        dg.edges[s,e]['cycle'] = G.edges[edge]['cycle']

                    if not G.edges[edge]['shadow']:
                        dg.edges[s,e]['weight'] = G.edges[edge]['weight']
                
                for raw in order.id_set: # only copy over actually used ids 
                    dg.nodes[s]['matches'][raw] = G.nodes[s]['matches'][raw]
                    dg.nodes[e]['matches'][raw] = G.nodes[e]['matches'][raw]
                    dg.edges[s,e]['matches'][raw] = G.edges[s,e]['matches'][raw]
                dg.edges[s,e]['orders'].append(order)
           
            if not 'consensus_direction' in cc.edges[edge]:
                continue
                    
            next_mod = direction_mod * cc.edges[edge]['consensus_direction'] # update modifier

            if other_exists:
                continue   # this is a revisit of an already used edge

            stack.append( (other, next_mod) )


############## Find and Assemble Contigs

####### Path Extraction

def topological_sort_reduction(G, O):

    found = set()

    indegree_map = {v: d for v, d in G.in_degree() if d > 0}
    zero_indegree = [v for v, d in G.in_degree() if d == 0]
    neighbour = set()

    if len(neighbour) == 0 and len(indegree_map) > 0:
        neighbour.add(next(iter(indegree_map)))

    while True:

        #print "Loop True ==========================="


        while zero_indegree:
            node = zero_indegree.pop()

            #print "Pop", node

            for _, child in G.out_edges(node):
                indegree_map[child] -= 1

                if indegree_map[child] == 0:
                    zero_indegree.append(child)
                    del indegree_map[child]
                    if child in neighbour:
                        neighbour.remove(child)
                else:
                    neighbour.add(child)
            #G.remove_node(node) DEBUG ONLY

        if not indegree_map:
            break # break global loop if no nodes are left


        #for edge in G.edges(): 
          #print edge, G.edges[edge]

        #print "Neighbour", neighbour
        #print "Indeg", len(indegree_map) 
        #print "ZeroDeg", len(zero_indegree) 
 
        empty_neighbour = False
        if len(neighbour) == 0:
            for s in indegree_map:
                neighbour.add(s)
            empty_neighbour = True

        min_score = -1;
        v_min = None
        for v in neighbour:
            score = 0
            for s, _ in G.in_edges(v):
                if s not in indegree_map:
                    continue
                score = score + 1 #G.edges[s, v_min]['weight']
            if v_min == None or score < min_score:
                min_score = score
                v_min = v

        for s, _ in G.in_edges(v_min):
            if s in indegree_map:
                G.remove_edge(s, v_min)
                O.edges[s, v_min]['shadow'] = True

        del indegree_map[v_min]
        zero_indegree.append(v_min)
        neighbour.remove(v_min)

        if empty_neighbour:
            neighbour.clear()

def find_cluster_weight(dg):

    #start_ = timer()

    topo = list(nx.topological_sort(dg))
    n_to_i = {}
    i_to_n = {}

    for i, n in enumerate(topo):    
        n_to_i[n] = i
        i_to_n[i] = n
        #print n, " to ", i

    for e in dg.edges():
        dg.edges[e]['cluster_weight'] = 0

    for nt in topo:

        dg.nodes[nt]['out'] = set()
        for e in dg.out_edges(nt):
            dg.nodes[nt]['out'].add(n_to_i[e[1]])

        dg.nodes[nt]['in'] = set()
        for e in dg.in_edges(nt):
            dg.nodes[nt]['in'].add(n_to_i[e[0]])

    for nt in topo:
  
        candidates = [ (dg.nodes[nt]['out'], [n_to_i[nt]]) ]  # initialize with full set 
        for id_ in sorted(dg.nodes[nt]['out']): # move through full list of out_going

            active = i_to_n[id_] 

            for in_ in dg.nodes[active]['in']:

                for open_, visited in candidates:
                    if visited[-1] == in_ and id_ in open_:
                        candidates.append( ( open_.intersection(dg.nodes[active]['out']), visited + [id_] ) )

            filtered = []
            for i in range(len(candidates)):
                dominated = False
                for j in range(len(candidates)):
                    if i != j and candidates[i][0].issubset(candidates[j][0]) and set(candidates[i][1]).issubset(set(candidates[j][1])):
                        dominated = True
                        break
                if not dominated:
                    filtered.append(candidates[i])
            candidates = filtered

        max_visited = []
        max_len = 0
        for  open_, visited in candidates:
            if len(visited) > max_len:
                max_visited = [visited]
                max_len = len(visited)
            elif len(visited) == max_len:
                max_visited.append(visited)

        for mv in max_visited:
            i = len(mv) - 1
            for u, v in zip(mv[:-1], mv[1:]):    
                nu = i_to_n[u]
                nv = i_to_n[v]
                dg.edges[nu,nv]['cluster_weight'] = dg.edges[nu,nv]['cluster_weight'] + i
                i = i - 1

    #end_ = timer()
    #print "Cluster", (end_ - start_)


def find_conservation_path(dg):
    #start_ = timer()

    topo = list(nx.topological_sort(dg))
    n_to_i = {}
    i_to_n = {}

    for i, n in enumerate(topo):    
        n_to_i[n] = i
        i_to_n[i] = n
        #print n, " to ", i

    finalized = []
    open_paths = [] 
    for n in topo:

        #print "N", n, n_to_i[n]

        if dg.out_degree(n) == 0: # True Endnode     
            for val, p in open_paths:
                if p[-1] == n:
                    finalized.append( [ x for x in p] )
                    finalized = [max(finalized,key=len)]
            continue

        max_outs = []
        max_out = 0
        for e in dg.out_edges(n):
            if dg.edges[e]['cluster_weight'] > max_out:
                max_out = dg.edges[e]['cluster_weight']
                max_outs = [ e ]
            elif dg.edges[e]['cluster_weight'] == max_out:
                max_outs.append(e)

        #print "Mo", max_outs, max_out

        max_ins = []
        max_in = 0
        for i, (val, p) in enumerate(open_paths):
            if p[-1] == n:
                if val > max_in:
                    max_in = val
                    max_ins = [ i ]
                elif val == max_in:
                    max_ins.append(i)

        #print "Mi", max_ins, max_in
        
        for mo in max_outs:
            if len(max_ins) > 0:
                #for mi in max_ins:
                    mi = max( max_ins, key=lambda x:len(open_paths[x][1]) )
                    open_paths.append( (max_out , open_paths[mi][1] + [mo[1]]) )
            else:
                open_paths.append( (max_out, [ mo[0], mo[1] ])  )

        
        open_paths[:] = [x for x in open_paths if n_to_i[x[1][-1]] > n_to_i[n]]
        #print "OP", len(open_paths) #" MEM ", process.memory_info().rss

    #end_ = timer()
    #print "Find Path", (end_ - start_)

    return finalized[0]

def iter_longest_path_extraction(dg):

    # initial copy
    dgc = dg.copy()
    
    # showdows should be ignored
    delete = []
    for e in dgc.edges():
        if dgc.edges[e]['shadow']:
            delete.append(e)
    for e in delete:
        dgc.remove_edge(e[0], e[1])

    #start_ = timer()
    topological_sort_reduction(dgc, dg)
    #end_ = timer()
    #print "Cycle", (end_ - start_)

    find_cluster_weight(dgc)

    # iteratitvely build up the path, longest to smallest
    path_set = []
    #vi = 1
    while len(dgc.edges()) > 0: # while we have edges
        
        #print "<<<<<<< Find Conservation <<<<<<<"
        #print "Find Conservation", process.memory_info().rss
        #print "Graph_Size", len(dgc.nodes()), len(dgc.edges())
        #for edge in dgc.edges(): 
        #    print edge, dgc.edges[edge]['cluster_weight']

        lp = find_conservation_path(dgc)
        #print lp

        if len(lp) < 10:
            in_vi = False
            for s , _ in dg.in_edges(lp[0]):
                if 'visit_id' in dg.nodes[s]:
                    in_vi = True
            out_vi = False
            for _ , t in dg.out_edges(lp[-1]):
                if 'visit_id' in dg.nodes[t]:
                    out_vi = True

            if (not in_vi and not out_vi) or ( (in_vi or out_vi) and len(lp) > 5 ):
                path_set.append(lp)
        else:
            path_set.append(lp)

        for v in lp:
            dg.nodes[v]['visit_id'] = True
            dgc.remove_node(v)

    for n in dgc.nodes():
        path_set.append( [n] )

    #print "Out", path_set
    return path_set

def rec_cc(colour_correction, cid):
    while colour_correction[cid] != cid:
        cid = colour_correction[cid]
    return cid

def linearize_directed_graph(dg):
    path_set = iter_longest_path_extraction(dg) # all paths to build a node cover
    # try if we can improve this by inversing bi-edges
    colour_correction = {}
    colour_to_length = {}
    for i, path in enumerate(path_set):
        for pos, node in enumerate(path):
            dg.nodes[node]["path_id"] = i
        colour_correction[i] = i
        colour_to_length[i] = len(path)

    potential_joins = []
    for e in dg.edges():
        if dg.edges[e]['shadow']:
            st = e[0]
            en = e[1]

            if 'path_id' not in dg.nodes[st] or 'path_id' not in dg.nodes[en]:
                continue

            c1 = dg.nodes[st]["path_id"]
            c2 = dg.nodes[en]["path_id"]

            l1_s = path_set[c1].index(st)
            l1_e = colour_to_length[c1] - l1_s - 1
            l2_s = path_set[c2].index(en)
            l2_e = colour_to_length[c2] - l2_s - 1

            if c1 != c2 and l1_e < l1_s and l2_s < l2_e: # remaining parts are larger
                potential_joins.append( ( l1_e + l2_s, e) )
    potential_joins.sort()

    shadow_joins = 0
    for dist, e in potential_joins:

        if dist > 3:
            break         

        #print "JOIN STEP", path_set

        st = e[0]
        en = e[1]
        i1 = dg.nodes[st]["path_id"]
        i2 = dg.nodes[en]["path_id"]
        c1 = rec_cc(colour_correction, i1)
        c2 = rec_cc(colour_correction, i2)
        if c1 == c2:
            continue

        try:
            l1 = path_set[c1].index(st)
            l2 = path_set[c2].index(en)
        except ValueError:
            # one was not found, so an overruled combine!
            continue

        #print st, en, "i:" , i1, i2, "c:", c1, c2, "l:", l1, l2

        # recalc dist to safeguard!
        l1_e = colour_to_length[c1] - l1 - 1
        l2_s = l2

        if l1_e + l2_s != dist:
            continue

        # join both paths
        shadow_joins = shadow_joins + 1
        path_set[c1] = path_set[c1][:l1+1] + path_set[c2][l2:]
        path_set[c2] = []

        colour_correction[c2] = colour_correction[c1] # make those the same
        colour_to_length[c1] = len(path_set[c1])
        colour_to_length[c2] = 0

    return [x for x in path_set if len(x) > 1], shadow_joins # short cut for now TODO  

####### Assembly and Consensus Calling

def cluster_anchors(dg, illuminaIdBase, edge_indices, edge_list, cluster_mod, id_to_overlap):

    ## build edge overlap to find maximal consistend cliques for anchor building
    ## directions are always consistent, just check overlaps
    aoog = nx.Graph()
    for ei_1 in edge_indices:
        aoog.add_node(ei_1)
        for ei_2 in edge_indices:

            if ei_1 == ei_2:
                break

            Q1 = dg.edges[edge_list[ei_1]]['matches'][illuminaIdBase].overlap
            Q2 = dg.edges[edge_list[ei_2]]['matches'][illuminaIdBase].overlap

            # the interval of the query both share
            overlapQ = (max(Q1[0], Q2[0]), min(Q1[1], Q2[1]))

            # and if the alignment coordinates within the query contig intersect:
            if overlapQ[0] <= overlapQ[1]:  # TODO: should we care for the size of the interval here ?  
                aoog.add_edge(ei_1, ei_2)

    ## optimal version
    #anchor_cliques = []
    #while not aoog.is_empty():
    #    opt_max = nxc_e.find_cliques(G)
    #    max_set = {}
    #    for om in opt_max:
    #        if len(om) > len(max_set):
    #            max_set = om
    #    anchor_cliques.append(om)
    #    for n in om:
    #        aoog.remove_node(n)

    ## heuristic version
    _, anchor_cliques = nxc_a.clique_removal(aoog)

    ## correct if not universal 
    for i, ac in enumerate(anchor_cliques):
             
        corrected_id = (illuminaIdBase, i)
        common_overlap = -1

        # build up the overlap
        for ei in ac: # ei = edge index

            cluster_mod[ei][illuminaIdBase] = i

            if common_overlap == -1:
                common_overlap = dg.edges[edge_list[ei]]['matches'][illuminaIdBase].overlap
            else:
                other = dg.edges[edge_list[ei]]['matches'][illuminaIdBase].overlap
                common_overlap = (max(common_overlap[0], other[0]), min(common_overlap[1], other[1]))

        id_to_overlap[corrected_id] = common_overlap

def compute_corrected_nano( match , overlap):

    nc_l = (overlap[0] -  match.illumina_range[0]) / match.ratio
    nc_r = (match.illumina_range[1] - overlap[1]) / match.ratio

    if match.direction < 0: # flip
        t = nc_l
        nc_l = nc_r
        nc_r = t

    #print "NC", nc_l, nc_r, " TO ", match.nano_range

    return ( match.nano_range[0] + nc_l, match.nano_range[1] - nc_r)


def get_illumina_sequence(illumina_id, illumina_overlap, direction):

    if direction < 0:
        return reverse_complement( find_sequence_illumina(illumina_id, illumina_overlap[0], illumina_overlap[1]) )
    else:
        return find_sequence_illumina(illumina_id, illumina_overlap[0], illumina_overlap[1])

def get_nano_sequence(nano_id, nano_region, direction):

    if direction < 0:
        return reverse_complement( find_sequence_nano(nano_id, nano_region[0], nano_region[1]) )
    else:
        return find_sequence_nano(nano_id, nano_region[0], nano_region[1])

def get_sequence_left_of_anchor(dg, nano, illumina_id, illumina_overlap, direction):

    if direction < 0:
        if dg.nodes[nano]['matches'][illumina_id].direction < 0: # single invert
            illu_seq = reverse_complement( find_sequence_illumina(illumina_id, dg.nodes[nano]['matches'][illumina_id].illumina_range[0], illumina_overlap[0] ) )
        else:
            illu_seq = find_sequence_illumina(illumina_id, illumina_overlap[1], dg.nodes[nano]['matches'][illumina_id].illumina_range[1] )

        nano_seq = find_sequence_nano( nano, G.nodes[nano]['matches'][illumina_id].nano_range[1], dg.nodes[nano]['nano_length']-1 )
        return reverse_complement( illu_seq + nano_seq)
    else:
        if dg.nodes[nano]['matches'][illumina_id].direction < 0: # single invert
            illu_seq = reverse_complement( find_sequence_illumina(illumina_id, illumina_overlap[1], dg.nodes[nano]['matches'][illumina_id].illumina_range[1]) )
        else:
            illu_seq = find_sequence_illumina(illumina_id, dg.nodes[nano]['matches'][illumina_id].illumina_range[0], illumina_overlap[0])

        nano_seq = find_sequence_nano( nano, 0, G.nodes[nano]['matches'][illumina_id].nano_range[0] )
        return nano_seq + illu_seq

def get_sequence_right_of_anchor(dg, nano, illumina_id, illumina_overlap, direction):

    if direction < 0:
        if dg.nodes[nano]['matches'][illumina_id].direction < 0: # single invert
            illu_seq = reverse_complement( find_sequence_illumina(illumina_id, illumina_overlap[1], dg.nodes[nano]['matches'][illumina_id].illumina_range[1] ) )
        else:
            illu_seq = find_sequence_illumina(illumina_id, dg.nodes[nano]['matches'][illumina_id].illumina_range[0], illumina_overlap[0] )

        nano_seq = find_sequence_nano( nano, 0, G.nodes[nano]['matches'][illumina_id].nano_range[0] )
        return reverse_complement( nano_seq + illu_seq )
    else:
        if dg.nodes[nano]['matches'][illumina_id].direction < 0: # single invert
            illu_seq = reverse_complement( find_sequence_illumina(illumina_id, dg.nodes[nano]['matches'][illumina_id].illumina_range[0], illumina_overlap[0]) )
        else:
            illu_seq = find_sequence_illumina(illumina_id, illumina_overlap[1], dg.nodes[nano]['matches'][illumina_id].illumina_range[1])

        nano_seq = find_sequence_nano( nano, G.nodes[nano]['matches'][illumina_id].nano_range[1], dg.nodes[nano]['nano_length']-1 )
        return illu_seq + nano_seq

def get_anchor_sequence(dg, nano, illumina_id, illumina_overlap, direction):

    d = direction * dg.nodes[nano]['matches'][illumina_id].direction

    if d < 0:
        return reverse_complement( find_sequence_illumina(illumina_id, illumina_overlap[0], illumina_overlap[1]) )
    else:
        return find_sequence_illumina(illumina_id, illumina_overlap[0], illumina_overlap[1])

def get_sequence_between_anchors(dg, nano, left_id, right_id, left_overlap, right_overlap, direction):
    # overlap -> illumina anchor region

    illu_r = dg.nodes[nano]['matches'][right_id].illumina_range
    illu_l = dg.nodes[nano]['matches'][left_id].illumina_range

    ratio_r = dg.nodes[nano]['matches'][right_id].ratio
    ratio_l = dg.nodes[nano]['matches'][left_id].ratio

    nano_r = dg.nodes[nano]['matches'][right_id].nano_range
    nano_l = dg.nodes[nano]['matches'][left_id].nano_range

    #print "IL", illu_l 
    #print "IR", illu_r 
    #print "NL", nano_l 
    #print "NR", nano_r 
    #print "DIR", direction 

    left_correction = 0
    right_correction = 0

    if direction < 0: # inversed nano
        left_correction = 0
        right_correction = 0

        offset_error = nano_r[1] - nano_l[0]
        if offset_error > 0:

            cn_l = compute_corrected_nano( dg.nodes[nano]['matches'][left_id], left_overlap )
            cn_r = compute_corrected_nano( dg.nodes[nano]['matches'][right_id], right_overlap )

            #print "OFFSET", cn_l, cn_r

            if cn_l[0] < cn_r[1]:  # the actual consensus overlaps overlap, get minus distance between both with no sequence
                return int(cn_l[0] - cn_r[1]), None

            # no nano part, shorten the illuminas greedy    
            if dg.nodes[nano]['matches'][left_id].direction < 0:
                available_left = (illu_l[1] - left_overlap[1]) / ratio_l
                left_correction = illu_l[1] - left_overlap[1]
            else:
                available_left = (left_overlap[0] - illu_l[0]) / ratio_l
                left_correction = left_overlap[0] - illu_l[0]

            if available_left > offset_error: # we can get all off in this one, so do it
                left_correction = int(offset_error * ratio_l)
                offset_error = 0
            else:
                offset_error -= available_left

            if dg.nodes[nano]['matches'][right_id].direction < 0:
                available_right = (right_overlap[0] - illu_r[0]) / ratio_r
                right_correction = right_overlap[0] - illu_r[0]    
            else:
                available_right = (illu_r[1] - right_overlap[1]) / ratio_r
                right_correction = illu_r[1] - right_overlap[1]

            if available_right > offset_error: # we can get the rest off in this one, so do it (should always be the case)
                right_correction = int(offset_error * ratio_r)


        if dg.nodes[nano]['matches'][right_id].direction < 0: # double invert
            illumina_left_seq = reverse_complement( find_sequence_illumina(right_id, illu_r[0] + right_correction, right_overlap[0]) )
        else:
            illumina_left_seq = find_sequence_illumina( right_id, right_overlap[1] , illu_r[1] - right_correction )

        if dg.nodes[nano]['matches'][left_id].direction < 0: # double invert
            illumina_right_seq = reverse_complement( find_sequence_illumina(left_id, left_overlap[1] , illu_l[1] - left_correction) )
        else:
            illumina_right_seq = find_sequence_illumina( left_id, illu_l[0] + left_correction, left_overlap[0] )

        nano_seq = find_sequence_nano( nano, nano_r[1], nano_l[0] )

        seq = reverse_complement( illumina_left_seq + nano_seq + illumina_right_seq )
        return len(seq), seq
    
    else:

        offset_error = nano_l[1] - nano_r[0]   ## nano part is empty on sequence extraction
        if offset_error > 0:

            cn_l = compute_corrected_nano( dg.nodes[nano]['matches'][left_id], left_overlap )
            cn_r = compute_corrected_nano( dg.nodes[nano]['matches'][right_id], right_overlap )

            #print "OFFSET", cn_l, cn_r

            if cn_l[1] > cn_r[0]:  # the actual consensus overlaps overlap, get minus distance between both with no sequence
                return int(cn_r[0] - cn_l[1]), None   # inversions also caugth here, but likely rare

            # no nano part, shorten the illuminas greedy    
            if dg.nodes[nano]['matches'][left_id].direction < 0:
                available_left = (left_overlap[0] - illu_l[0]) / ratio_l
                left_correction = left_overlap[0] - illu_l[0]
            else:
                available_left = (illu_l[1] - left_overlap[1]) / ratio_l
                left_correction = illu_l[1] - left_overlap[1]

            if available_left > offset_error: # we can get all off in this one, so do it
                left_correction = int(offset_error * ratio_l)
                offset_error = 0
            else:
                offset_error -= available_left

            if dg.nodes[nano]['matches'][right_id].direction < 0:
                available_right = (illu_r[1] - right_overlap[1]) / ratio_r
                right_correction = illu_r[1] - right_overlap[1]
            else:
                available_right = (right_overlap[0] - illu_r[0]) / ratio_r
                right_correction = right_overlap[0] - illu_r[0]

            if available_right > offset_error: # we can get the rest off in this one, so do it (should always be the case)
                right_correction = int(offset_error * ratio_r)

        if dg.nodes[nano]['matches'][left_id].direction < 0: # single invert
            illumina_left_seq = reverse_complement( find_sequence_illumina(left_id, illu_l[0] + left_correction, left_overlap[0]) )
        else:
            illumina_left_seq = find_sequence_illumina( left_id, left_overlap[1] , illu_l[1] - left_correction )
        
        if dg.nodes[nano]['matches'][right_id].direction < 0: # single invert
            illumina_right_seq = reverse_complement( find_sequence_illumina(right_id, right_overlap[1] , illu_r[1] - right_correction) )
        else:             
            illumina_right_seq = find_sequence_illumina( right_id, illu_r[0] + right_correction, right_overlap[0] )

        nano_seq = find_sequence_nano( nano, nano_l[1], nano_r[0] )

        seq = illumina_left_seq + nano_seq + illumina_right_seq
        return len(seq), seq


def align_anchor_region(dg, ag, left_anchor, right_anchor, left_overlap, right_overlap):

    # overlap -> illumina anchor region with
    # illumina match id with a special index for disconnects!
    sequences = [] # array of sequences to collect, align and call consensus

    ## get all sequences
    g_dist = None
    for nano in ag.edges[left_anchor, right_anchor]['nanopores']:   
        direction = dg.nodes[nano]['direction']

        dist, seq = get_sequence_between_anchors(dg, nano, left_anchor[0][0], right_anchor[0][0], left_overlap, right_overlap, direction)  

        if seq != None:
            sequences.append(seq)
        if g_dist == None: ## take first distance as goal, so we are guaranteed a base sequence
            g_dist = dist

    ## add in as region
    ag.edges[left_anchor, right_anchor]['sequences'] = sequences
    ag.edges[left_anchor, right_anchor]['distance'] = g_dist

def update_consensus_base( old_seq, old_p1, old_p2, add_seq, add_p1, add_p2):

    if old_seq == None:
        return (add_seq, add_p1, add_p2)

    new_seq = old_seq
    if add_p1 < old_p1:
        new_seq = add_seq[: old_p1 - add_p1] + new_seq
    if add_p2 > old_p2:
        new_seq = new_seq + add_seq[(-1*(add_p2 - old_p2)):]

    return ( new_seq, min(old_p1, add_p1), max(old_p2, add_p2) )

def ordered_visit(adg, order, tap, id_to_overlap, start_node):

    sequence = None
    left_border = 0
    right_border = 0

    edge_queue = SortedSet( key=lambda x:(x[0],-x[1]) )
    node_queue = SortedSet()


    node_queue.add(adg.nodes[start_node]['order'])
    while len(node_queue) > 0:
        i = node_queue[0]
        v = order[i]
        node_queue.pop(0)

        if not 'visited' in adg.nodes[v]:

            adg.nodes[v]['visited'] = True
            for edge in adg.out_edges(v):
                edge_queue.add( (adg.nodes[edge[1]]['order'] , i) )
                node_queue.add( adg.nodes[edge[1]]['order'] )

            while len(edge_queue) > 0 and edge_queue[0][0] == i: # sweep up all edges ending in current one
                
                left_anchor = order[edge_queue[0][1]]
                right_anchor = order[edge_queue[0][0]]
 
                has_l = left_anchor in tap
                has_r = right_anchor in tap

                ovl_1 = id_to_overlap[left_anchor[0]]
                ovl_2 = id_to_overlap[right_anchor[0]]

                offset = adg.edges[left_anchor, right_anchor]['distance']

                l_len = ovl_1[1] - ovl_1[0] + 1
                r_len = ovl_2[1] - ovl_2[0] + 1

                if has_l and not has_r:

                    rpos = tap[left_anchor][1]
                    tap[right_anchor] = (rpos + offset + 1 , rpos + offset + 1 + r_len -1)

                    if offset > 0 :
                        sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.edges[left_anchor, right_anchor]['sequences'][0], rpos + 1, rpos + offset)
                    sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.nodes[right_anchor]['anchor_seq'] , tap[right_anchor][0], tap[right_anchor][1])

                elif not has_l and has_r:

                    rpos = tap[right_anchor][0]
                    tap[left_anchor] = (rpos - offset - 1 - (l_len - 1), rpos - offset - 1 )   

                    if offset > 0 :
                        sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.edges[left_anchor, right_anchor]['sequences'][0], rpos - offset, rpos)
                    sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.nodes[left_anchor]['anchor_seq'] , tap[left_anchor][0], tap[left_anchor][1])

                elif not has_l and not has_r:

                    tap[left_anchor] = (0, l_len - 1)
                    tap[right_anchor] = (l_len + offset, l_len + offset + r_len - 1 )

                    if offset > 0:
                        sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.edges[left_anchor, right_anchor]['sequences'][0], l_len, l_len + offset - 1)

                    sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.nodes[left_anchor]['anchor_seq'] , tap[left_anchor][0], tap[left_anchor][1])
                    sequence, left_border, right_border = update_consensus_base(sequence, left_border, right_border, adg.nodes[right_anchor]['anchor_seq'] , tap[right_anchor][0], tap[right_anchor][1])

                # true true => longer edge we do not care for
                 
                edge_queue.pop(0)

        else:
            while len(edge_queue) > 0 and edge_queue[0][0] == i:
                edge_queue.pop(0)

    return (sequence, left_border, right_border)


def assemble_path(path, dg, out_id, of_query, of_paf, of_target, of_log):

    #print "Assemble", out_id
    of_log.write("Assemble Path 'Prokrastinator_"+ out_id + "' of length " + str(len(path)) + ".\n")
    #print "Assemble Path 'Prokrastinator_"+ out_id + "' of length " + str(len(path))
    #print path

    Candidate = collections.namedtuple('Candidate', 'open_ids visited_ids score kinks edge_list order_list modifier')
    # max bit, min kink

    #start_ = timer()
    ####### Build up the Anchor Connection DAG
    candidate_list = []
    candidate_list.append( Candidate( open_ids = set(), visited_ids = set(), score = 0, kinks = 0, edge_list = [], order_list = [], modifier = [] )  ) 

    for v, u in zip(path[:-1], path[1:]): 

        e = (v,u)

        #print "Path", e        

        next_candidate_list = []
        for order in dg.edges[e]['orders']:

            sub_candidates = []
            for candidate in candidate_list:

                ## core is easy
                bs = candidate.score + order.score

                ## retrieve the edge
                id_set = order.id_set
                base_node = order.base_node

                if dg.nodes[base_node]['direction'] < 0: # inversed base
                    id_set.reverse()

                ## loop and check all anchors
                new_kinks = candidate.kinks
                next_open_ids = set()
                edge_modifier = set()
                for a in id_set: # these are in order, we CANNOT have inversions to open ids in theory...
                    if a not in candidate.open_ids and a in candidate.visited_ids:
                        edge_modifier.add(a)
                        new_kinks += 1
                    next_open_ids.add(a)

                sub_candidates.append( Candidate( open_ids = next_open_ids, visited_ids = candidate.visited_ids | set(id_set) , score = bs, kinks = new_kinks, edge_list = candidate.edge_list + [e], order_list = candidate.order_list + [order], modifier = candidate.modifier + [edge_modifier] ) )

            ## choice function
            min_kin = -1
            max_score = 0
            for candidate in sub_candidates:
                if min_kin == -1 or candidate.kinks < min_kin or ( candidate.kinks == min_kin and candidate.score > max_score ):
                    min_kin = candidate.kinks
                    max_score = candidate.score
            next_candidate_list = next_candidate_list + [ x for x in sub_candidates if x.kinks == min_kin and x.score == max_score ]
        candidate_list = next_candidate_list

    ## final best
    min_kin = -1
    max_score = 0
    for candidate in candidate_list:
        if min_kin == -1 or candidate.kinks < min_kin or ( candidate.kinks == min_kin and candidate.score > max_score ):
            min_kin = candidate.kinks
            max_score = candidate.score
    best_candidate =  [x for x in candidate_list if x.kinks == min_kin and x.score == max_score ][0] # choose one at random
     
    #end_ = timer()
    #print "Sub Candidates-Path", (end_ - start_)

    #start_ = timer()
    # cluster each anchor row!
    cluster_collector = {}
    cluster_mod = [ {} for i in range(len(best_candidate.edge_list)) ]
    for i in range(len(best_candidate.edge_list)):

        order = best_candidate.order_list[i]
        edge = best_candidate.edge_list[i]
        id_set = order.id_set

        #for ci in cluster_collector:        
        #    if ci not in id_set: # end of strech found
        #        # call cluster
        #        cluster_anchors(dg, ci, cluster_collector[ci], best_candidate.edge_list, cluster_mod, id_to_overlap)
                
        for match in id_set:
            cluster_collector.setdefault(match, []).append(i)

    for ci in cluster_collector:
        cluster_anchors(dg, ci, cluster_collector[ci], best_candidate.edge_list, cluster_mod, id_to_overlap)

    #end_ = timer()
    #print "Sub Clustering", (end_ - start_)


    #start_ = timer()
    # add info to node list
    nodes_info = [ [] for i in range(len(best_candidate.edge_list) + 1) ]
    nodes = [ "" for i in range(len(best_candidate.edge_list) + 1) ]
    match_modifier = {}
    for i in range(len(best_candidate.edge_list)):

        ## base info
        order = best_candidate.order_list[i]
        edge = best_candidate.edge_list[i]
        modifier = best_candidate.modifier[i]

        ## update modifier
        for mod in modifier:
            match_modifier[mod] = match_modifier.setdefault(mod, 0) + 1

        ## ID set
        id_set = order.id_set
        base_node = order.base_node
        if dg.nodes[base_node]['direction'] < 0: # inversed base
            id_set.reverse()

        for match in id_set:
            match_co = ((match, cluster_mod[i][match]), 0)
            if match in match_modifier:
                match_co = ((match, cluster_mod[i][match]), match_modifier[match])

            nano_a = dg.nodes[edge[0]]['matches'][match].nano_range
            nodes_info[i].append( (nano_a, match_co) )

            nano_b = dg.nodes[edge[1]]['matches'][match].nano_range
            nodes_info[i + 1].append( (nano_b, match_co) )
           
        ## Nodes
        nodes[i] = edge[0] # n-1 overwrite, but faster than if
        nodes[i+1] = edge[1]

    pre_log = {}
    post_log = {}

    adg = nx.DiGraph()
    for i in range(len(nodes)):
          
        node = nodes[i]

        def cmp_ni(a, b):
            if a[0] == b[0]: # same region/same nano!
                if dg.nodes[node]['matches'][a[1][0][0]].direction < 0:
                    return cmp(id_to_overlap[b[1][0]] , id_to_overlap[a[1][0]])
                else:
                    return cmp(id_to_overlap[a[1][0]] , id_to_overlap[b[1][0]])
            return cmp(a[0], b[0])

        nodes_info[i].sort(cmp_ni) # sort by range

        if dg.nodes[node]['direction'] < 0: # inversed base
            nodes_info[i].reverse()

        of_log.write("Node "+ node + " Direction " + str(dg.nodes[node]['direction']) + " " + str(nodes_info[i]) + ".\n")
        #print "Node "+ node + " Direction " + str(dg.nodes[node]['direction']) + " " + str(nodes_info[i]) 

        ni = nodes_info[i]   

        if len(ni) == 0:
            continue

        last = ni[0][1]
        last_nr = ni[0][0]
        for nr, match in ni:

            if not adg.has_node(match):
                adg.add_node(match)
                adg.nodes[match]['anchor_seq'] = get_anchor_sequence(dg, node, match[0][0], id_to_overlap[match[0]], dg.nodes[node]['direction'])

            if match == last:
                continue

            if not adg.has_node(last):
                adg.add_node(last)
                adg.nodes[last]['anchor_seq'] = get_anchor_sequence(dg, node, last[0][0], id_to_overlap[last[0]], dg.nodes[node]['direction'])

            edge = (last, match)
            if ( last_nr[1] > nr[1] and last_nr[0] < nr[0] ) or ( last_nr[1] < nr[1] and last_nr[0] > nr[0] ): # we have a case of containment!

                cn_l = compute_corrected_nano( dg.nodes[node]['matches'][last[0][0]], id_to_overlap[last[0]] )
                cn_r = compute_corrected_nano( dg.nodes[node]['matches'][match[0][0]], id_to_overlap[match[0]] )

                if ( dg.nodes[node]['direction'] > 0  and ( cn_l[0] > cn_r[0] or (cn_l[0] == cn_r[0] and cn_l[1] > cn_r[1])) ) \
                   or ( dg.nodes[node]['direction'] < 0  and ( cn_l[0] < cn_r[0] or (cn_l[0] == cn_r[0] and cn_l[1] < cn_r[1])) ):
                    edge = (match, last)

            if not adg.has_edge(edge[0], edge[1]):
                adg.add_edge(edge[0], edge[1])

            #print "Add Edge", edge           

            adg.edges[ edge ].setdefault('nanopores', []).append( node ) 

            last = match
            last_nr = nr
        
        first_id = ni[0][1]
        second_id = ni[-1][1]

        adg.nodes[first_id].setdefault('pre_sequences', []).append( get_sequence_left_of_anchor(dg, node, first_id[0][0], id_to_overlap[first_id[0]], dg.nodes[node]['direction']) )
        adg.nodes[second_id].setdefault('post_sequences', []).append( get_sequence_right_of_anchor(dg, node, second_id[0][0], id_to_overlap[second_id[0]], dg.nodes[node]['direction']) ) 
    
        pre_log[node] = len(get_sequence_left_of_anchor(dg, node, first_id[0][0], id_to_overlap[first_id[0]], dg.nodes[node]['direction']))
        post_log[node] = len(get_sequence_right_of_anchor(dg, node, second_id[0][0], id_to_overlap[second_id[0]], dg.nodes[node]['direction']))

    ####### Find ordering of the DAG, and build up the anchor_position

    for edge in adg.edges():
        left_anchor = edge[0]
        right_anchor = edge[1]

        align_anchor_region(dg, adg, left_anchor, right_anchor, id_to_overlap[left_anchor[0]], id_to_overlap[right_anchor[0]])

    order = []
    order_g = nx.topological_sort(adg)
    for i, n in enumerate(order_g):
        adg.nodes[n]['order'] = i
        order.append(n)
    
    tap = {} # total_anchor_positions

    global_sequence, g_pos1, g_pos2 = ordered_visit(adg, order, tap, id_to_overlap, order[0])  # we create a base offset of all nodes reachable from order[0]          
    
    if len(adg.nodes()) == 1:
        anchor = list(adg.nodes())[0]
        ovl = id_to_overlap[anchor[0]]
        tap[anchor] = (0, ovl[1] - ovl[0])
        global_sequence =  adg.nodes[anchor]['anchor_seq']
        g_pos1 = 0 
        g_pos2 = ovl[1] - ovl[0]

    additional_paths = []
    path_added = []
    for node in order[1:]:  

        if 'visited' in adg.nodes[node]:
            continue

        tap2 = {} # total_anchor_positions 
        local_sequence, l_pos1, l_pos2 = ordered_visit(adg, order, tap2, id_to_overlap, node)

        if len(tap2) == 0:
            anchor = node
            ovl = id_to_overlap[anchor[0]]
            tap2[anchor] = (0, ovl[1] - ovl[0])
            local_sequence =  adg.nodes[anchor]['anchor_seq']
            l_pos1 = 0 
            l_pos2 = ovl[1] - ovl[0]

        additional_paths.append( (local_sequence, l_pos1, l_pos2, tap2) )
        path_added.append(False)

    loop = True
    while loop:
        loop = False

        for i, (local_sequence, l_pos1, l_pos2, tap2) in enumerate(additional_paths):

            if path_added[i] == True:
                continue

            for match in tap2: # find connection to existing "scaffolding"
                found = False
                for edge in adg.out_edges(match):
                    if edge[1] in tap:  # we found a connection
                        group_offset = tap[edge[1]][0] - adg.edges[edge]['distance'] - tap2[edge[0]][1] - 1

                        if len(adg.edges[edge]['sequences']) > 0:
                            local_sequence, l_pos1, l_pos2 = update_consensus_base(local_sequence, l_pos1, l_pos2, adg.edges[edge]['sequences'][0], tap2[edge[0]][1] + 1, tap2[edge[0]][1] + adg.edges[edge]['distance'])
                        found = True
                        break
                if found:
                    break
                for edge in adg.in_edges(match):
                    if edge[0] in tap:  # we found a connection
                        group_offset = tap[edge[0]][1] + adg.edges[edge]['distance'] + 1 - tap2[edge[1]][0]

                        if len(adg.edges[edge]['sequences']) > 0:
                            local_sequence, l_pos1, l_pos2 = update_consensus_base(local_sequence, l_pos1, l_pos2, adg.edges[edge]['sequences'][0], tap2[edge[1]][0] - adg.edges[edge]['distance'] , tap2[edge[1]][0] - 1)
                        found = True
                        break
                if found:
                    break

            if not found:
                loop = True

                continue

            path_added[i] = True
            for match in tap2: # copy over corrected to tap
                tap[match] = (tap2[match][0] + group_offset, tap2[match][1] + group_offset)
            global_sequence, g_pos1, g_pos2 = update_consensus_base(global_sequence, g_pos1, g_pos2, local_sequence, l_pos1 + group_offset, l_pos2 + group_offset)

    for node in adg.nodes():
        if 'pre_sequences' in adg.nodes[node]:
            max_seq = max( adg.nodes[node]['pre_sequences'], key=lambda x:len(x))
            global_sequence, g_pos1, g_pos2 = update_consensus_base(global_sequence, g_pos1, g_pos2, max_seq, tap[node][0] - len(max_seq) , tap[node][0] - 1)

        if 'post_sequences' in adg.nodes[node]:
            max_seq = max( adg.nodes[node]['post_sequences'], key=lambda x:len(x))
            global_sequence, g_pos1, g_pos2 = update_consensus_base(global_sequence, g_pos1, g_pos2, max_seq, tap[node][1] + 1, tap[node][1] + len(max_seq) )

    glmp = -1 * g_pos1 # global_left_most_postion , correction_offset for all output!
    
    of_log.write("TAP "+ str(tap) + ".\n")
    #print "TAP "+ str(tap)
    #print "glmp", glmp
    for node in adg.nodes():        
        of_log.write("ID_TO_OVRL " + str(node) + " : " + str(id_to_overlap[node[0]]) + ".\n")
        #print "ID_TO_OVRL " + str(node) + " : " + str(id_to_overlap[node[0]])

    for ni in range(len(nodes_info)):
        of_log.write( "REG " + str(nodes[ni])+ " : " + str(tap[nodes_info[ni][0][1]][0] - pre_log[nodes[ni]] + glmp) + " "+  str(tap[nodes_info[ni][-1][1]][1] + post_log[nodes[ni]] + glmp) + "\n")
        #print  "REG " + str(nodes[ni])+ " : " + str(tap[nodes_info[ni][0][1]][0] - pre_log[nodes[ni]] + glmp) + " "+  str(tap[nodes_info[ni][-1][1]][1] + post_log[nodes[ni]] + glmp) + " ; " + str(tap[nodes_info[ni][0][1]][0] - pre_log[nodes[ni]]) + " "+  str(tap[nodes_info[ni][-1][1]][1] + post_log[nodes[ni]])

    #if len(global_sequence) != g_pos2 - g_pos1 +1:
    #    print "WARNING"

    #end_ = timer()
    #print "Sub Region Extract", (end_ - start_)

    # target
    TN = "Prokrastinator_"+str(out_id)
    GSL = str(len(global_sequence))
    of_target.write(">"+TN+"\n")
    write_sequence(of_target, global_sequence)

    # query and PAF
    query_index = 0
    for edge in adg.edges(): 

        left_anchor = edge[0]
        right_anchor = edge[1]

        for seq in adg.edges[left_anchor, right_anchor]['sequences']:

            if len(seq) == 0:
                continue

            QSN = "Middle."+str(out_id)+"."+str(query_index)

            of_query.write(">"+QSN+"\n")
            write_sequence(of_query, seq)

            LB = tap[left_anchor][1] + 1 + glmp
            RB = tap[right_anchor][0] - 1 + glmp

            of_paf.write(QSN+"\t"+str(len(seq))+"\t0\t"+str(len(seq))+"\t+\t"+TN+"\t"+ GSL+ "\t" + str(LB) + "\t" + str(RB) + "\t"+ str(RB -LB +1)+ "\t"+ str(RB -LB +1) + "\t255" +"\n")

            query_index = query_index + 1

            of_log.write("Mid "+ str(LB) + "-" + str(RB) + " : " + str(left_anchor) + ", " + str(right_anchor) + " ; len=" + str(len(seq)) + ", stretch=" + str(RB-LB+1) + "\n")

    for node in adg.nodes():


        #of_query.write(">"+"QP_Anchor_"+str(out_id)+"_"+str(query_index)+"\n")
        #write_sequence(of_query, adg.nodes[node]['anchor_seq'])
        #query_index = query_index + 1

        if 'pre_sequences' in adg.nodes[node]: 
          for seq in adg.nodes[node]['pre_sequences']:

            if len(seq) < 200:
                continue

            QSN = "Left."+str(out_id)+"."+str(query_index)

            of_query.write(">"+QSN+"\n")
            write_sequence(of_query, seq)

            RB = tap[node][0] - 1 + glmp
            LB = RB - len(seq) + 1

            of_paf.write(QSN+"\t"+str(len(seq))+"\t0\t"+str(len(seq))+"\t+\t"+TN+"\t"+ GSL+ "\t" + str(LB) + "\t" + str(RB) + "\t"+ str(RB -LB +1)+ "\t"+ str(RB -LB +1) + "\t255" +"\n")

            query_index = query_index + 1

            of_log.write("Pre "+ str(LB) + "-" + str(RB) + " : " + str(node) + " ; len=" + str(len(seq)) + ", stretch=" + str(RB-LB+1) + "\n")

        if 'post_sequences' in adg.nodes[node]:
          for seq in adg.nodes[node]['post_sequences']:

            if len(seq) < 200:
                continue

            QSN = "Right."+str(out_id)+"."+str(query_index)

            of_query.write(">"+QSN+"\n")
            write_sequence(of_query, seq)

            LB = tap[node][1] + 1 + glmp
            RB = LB + len(seq) - 1

            of_paf.write(QSN+"\t"+str(len(seq))+"\t0\t"+str(len(seq))+"\t+\t"+TN+"\t"+ GSL+ "\t" + str(LB) + "\t" + str(RB) + "\t"+ str(RB -LB +1)+ "\t"+ str(RB -LB +1) + "\t255" +"\n")

            query_index = query_index + 1

            of_log.write("Post "+ str(LB) + "-" + str(RB) + " : " + str(node) + " ; len=" + str(len(seq)) + ", stretch=" + str(RB-LB+1) + "\n")

    #### add in contained

    for i in range(len(nodes)):

        node = nodes[i]
        infos = nodes_info[i]

        info_map = {} # map original ID to local anchor
        for info in infos:
            info_map[ info[1][0][0] ] = info[1]

        #print "InfoMap: ", node, info_map

        for ce in dg.nodes[node]['contain_elements']:

            #print "Contain_Element: ", ce
            #print "In: ", node, infos

            contain_info = []
            for mid in ce.matches:
                if mid in info_map:
                    contain_info.append( (ce.matches[mid].nano_range, mid) )

            if len(contain_info) == 0:
                continue
            
            contain_info.sort()
            direction = ce.direction * dg.nodes[node]['direction']
            if direction < 0:
               contain_info.reverse()

            #print "CI", contain_info, direction

            global_ranges = []
            for u in contain_info:
                tap_id = info_map[u[1]]
                tap_dir = dg.nodes[node]['matches'][u[1]].direction * dg.nodes[node]['direction']
                if tap_dir < 0:
                    illu_ref = id_to_overlap[tap_id[0]][0]
                else:
                    illu_ref = id_to_overlap[tap_id[0]][1]

                total_ref = tap[tap_id][1] + glmp

                #print "C", tap_id
                #print "ref", illu_ref
                #print "dirs CEC", ce.direction, "CEM", ce.matches[u[1]].direction, "OM", dg.nodes[node]['matches'][u[1]].direction, "N", dg.nodes[node]['direction']
                #print "mod_from", ce.matches[u[1]].illumina_range
                #print "Total_Tap", total_ref

                cont_dir = ce.matches[u[1]].direction * direction
                if cont_dir < 0:
                    offset = ce.matches[u[1]].illumina_range[0] - illu_ref
                    global_ranges.append( (total_ref - offset - (ce.matches[u[1]].illumina_range[1] - ce.matches[u[1]].illumina_range[0]), total_ref - offset) ) 
                else:
                    offset = ce.matches[u[1]].illumina_range[1] - illu_ref
                    global_ranges.append( (total_ref + offset - (ce.matches[u[1]].illumina_range[1] - ce.matches[u[1]].illumina_range[0]), total_ref + offset) )                   

                #print "R", cont_dir, global_ranges[-1]

            seqs = []
            for i in range(len(global_ranges)):

                gr = global_ranges[i]
                illu_id = contain_info[i][1]
                
                illu_range = ce.matches[illu_id].illumina_range
                nano_range = ce.matches[illu_id].nano_range
                illu_direction = ce.matches[illu_id].direction
               
                seq = get_illumina_sequence(illu_id, (illu_range[0], illu_range[1]), illu_direction * direction )
                seqs.append( (seq, gr[0], gr[1], "Illumina_Match") )


                if i == len(global_ranges) - 1:
                    if direction < 0:
                        seq = get_nano_sequence(ce.nano, (0, nano_range[0] - 1), direction)
                    else:
                        seq = get_nano_sequence(ce.nano, (nano_range[1], ce.nano_length - 1), direction)
                    #seqs.append( (seq, gr[1] + 1, gr[1] + len(seq), "Left") )
                if i == 0:
                    if direction < 0:
                        seq = get_nano_sequence(ce.nano, (nano_range[1], ce.nano_length - 1), direction)
                    else:
                        seq = get_nano_sequence(ce.nano, (0, nano_range[0] - 1), direction)
                    #seqs.append( (seq, gr[0] - len(seq), gr[0] - 1, "Right") )
                    continue

                pre_nano =  ce.matches[contain_info[i - 1][1]].nano_range
                seq = get_nano_sequence(ce.nano, (pre_nano[1] + 1, nano_range[0] - 1), direction)
                seqs.append( (seq, global_ranges[i-1][1] + 1, gr[0] - 1, "Nano_Middle") )

            for s in seqs:  

                if len(s[0]) < 200:
                    continue

                QSN = "Contain_" + s[3] +"." +str(out_id)+"."+str(query_index)

                of_query.write(">"+QSN+"\n")
                write_sequence(of_query, s[0])

                LB = s[1]
                RB = s[2]

                of_paf.write(QSN+"\t"+str(len(s[0]))+"\t0\t"+str(len(s[0]))+"\t+\t"+TN+"\t"+ GSL+ "\t" + str(LB) + "\t" + str(RB) + "\t"+ str(RB -LB +1)+ "\t"+ str(RB -LB +1) + "\t255" +"\n")
                query_index = query_index + 1

                of_log.write("Contain "+ str(LB) + "-" + str(RB) + " " + node + "\n")


############## Plotting for debug

def plot_graph(dg, i, paths):

    pos = nx.circular_layout(dg)

    scale = 3
    maxi = 0
    clist = ['red', 'green', 'blue', 'yellow', 'magenta', 'orange', 'cyan', 'purple']
    for ci, path in enumerate(paths):
        for u, v in zip(path[:-1], path[1:]):
             dg.edges[u,v]['color'] = clist[ci]

        slice_c = math.pi / float(len(path))
        for ui, u in enumerate(path):
            angle = slice_c * ui;
            newX = 20 + 10 * math.cos(angle)
            newY = 2*(ci + 1) + 10 * math.sin(angle)

            dg.nodes[u]['pos'] = [newX*scale, newY*scale]
 
    for e in dg.edges():
        dg.edges[e]['label'] = 'S' if dg.edges[e]['shadow'] else ''

    todo = []
    for n in dg.nodes():
        if 'pos' not in dg.nodes[n]:
            todo.append(n)
    slice_c = math.pi / float(max(len(todo), 1))
    for ui, n in enumerate(todo):
        angle = slice_c * ui;
        newX = 20 + 10 * math.cos(angle)
        newY = 0 + 10 * math.sin(angle) 
           
        dg.nodes[n]['pos'] = [newX*scale, newY*scale]

    colors = [dg.edges[e]['color'] if 'color' in dg.edges[e] else 'black' for e in dg.edges()]
    pos = nx.get_node_attributes(dg,'pos')
    labels = nx.get_edge_attributes(dg,'label')

    nx.draw(dg, pos, node_color='black', edge_color= colors, node_size=210,
                    with_labels=True, font_color="black", font_size=6)
    nx.draw_networkx_edge_labels(dg, pos, edge_labels=labels)

    plt.savefig(tmp_dir + "/prokrast_" + str(i) + ".png")
    plt.close() 


############## Core Loop

def custom_plain_bfs(G, source):
    seen = set()
    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                yield v
                seen.add(v)
                for edge in G.edges(v):
                    if 'consensus_direction' in G.edges[edge]:
                        nextlevel.add(edge[1])

def custom_connected_components(G):
    seen = set()
    for v in G:
        if v not in seen:
            c = set(custom_plain_bfs(G, v))
            yield c
            seen.update(c)



total = len(list(custom_connected_components(G)))

assembly_index = 0
for i, ccc in enumerate(custom_connected_components(G)):

    cc = G.subgraph(ccc)

    max_value = 0
    max_node = -1
    for node in cc.nodes(): # find longest edge
        if cc.nodes[node]['nano_length'] > max_value:
            max_value = cc.nodes[node]['nano_length']
            max_node = node

    #print "UNDIRECTED =================================="
    #for e in cc.edges():
    #    print e, cc.edges[e]['shadow']

    print "------------------- Graph "+str(i+1)+"/" + str(total)

    #start = timer()

    dg = nx.DiGraph() # create new directed graph, we start empty
    direction_graph(cc, max_node, dg)

    #end = timer()
    #print "Directing Graph", (end - start)

   
    #print "DIRECTED >>>>>>>>>>>>>>>>>>>>>>"
    #for e in dg.edges():
    #    print e, dg.edges[e]['orders'], dg.edges[e]['shadow']
    #print "DIRECTED <<<<<<<<<<<<<<<<<<<<<<"

    of_log.write("Linearize Graph of " + str(i) + " of size " + str(len(dg.nodes()))) 
    #start = timer()
    paths, shadow_joins = linearize_directed_graph(dg)
    #end = timer()
    #print "Linearize", (end - start)

    of_log.write("Graph "+str(i)+"/" + str(total) + " of size " + str(len(dg.nodes())) + " seperated into " + str(len(paths)) + " paths. \n")
    of_log.write("Shadow-Joins " + str(shadow_joins) + " \n")

    #print "Linearized", len(paths)

    #if len(paths) < 9:
    #    plot_graph(dg, i, paths)

    for path in paths:
        #start = timer()
        assemble_path(path, dg, str(assembly_index), of_query, of_paf, of_target, of_log)
        #end = timer()
        #print "Assemble", (end - start)

        assembly_index = assembly_index + 1

of_query.close()
of_paf.close()
of_target.close()
of_log.close()


