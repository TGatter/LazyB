#!/usr/bin/env python

import networkx as nx
import sys, os
from Bio import SeqIO
import random
import time

start_time = time.time()

paf_anchors = sys.argv[1] # anchors to nanopores
seq_file_base = sys.argv[2] # nanopore reads
output_file = sys.argv[3] # scrubbed nanopore reads
tmp_dir = sys.argv[4]
output = open(output_file, "w")

subset_size = 60000

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

# set minimum length for mapping hits
threshold_length = 500
# create graph
G = nx.Graph()

prevHitId = ""  # just to prevent the error at the first line
chunkNodes = {}
############## Read Blast in and build base graph

for line in open(paf_anchors, "r"):

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
#11 	int 	Mapping quality (0-255; 255 for missing)

    linetemp = line.rstrip().split("\t")
    if len(linetemp) == 1:
        continue

    id_1, id_2 = linetemp[0], linetemp[5]  # remember about 0-based counting

    len_1, len_2 = int(linetemp[1]), int(linetemp[6])  # remember about 0-based counting

    s_1, s_2 = int(linetemp[2]), int(linetemp[7])
    e_1, e_2 = int(linetemp[3]), int(linetemp[8])

    if e_1 - s_1 < 500:         # only hits >= 500 bp count
        continue

    
    # add node to the graph if necessary 
    if not G.has_node(id_2):  # if this node does not exist, then create it
        G.add_node(id_2)#, print_index=node_index)
        G.nodes[id_2]['length'] = len_2 # set length of nanopore read 
        G.nodes[id_2]['illu_to_ranges'] = {} # initialize empty matches to anchors
        G.nodes[id_2]['seq_to_ranges'] = {} # initialize empty matches to other nanos

    # add the anchor hit range to the vertex information  
    # we exclude information of second hits of the same anchor-nanopore pair
    if id_1 in G.nodes[id_2]['illu_to_ranges']:
        continue  

    G.nodes[id_2]['illu_to_ranges'][id_1] = (s_2, e_2)
    
    # are we still in the same scaffold? If not, let's start a new one
    if id_1 != prevHitId:
        chunkNodes = []
        prevHitId = id_1

    # here are reads associated with the same anchor
    # and now iterate only through these reads
    for prevId in chunkNodes:
        # add edge if necessary
        if not G.has_edge(prevId, id_2):
            G.add_edge(prevId, id_2)
    chunkNodes.append(id_2)  

total_nodes = len(G.nodes())
print("\n" + "# Nodes " + str(total_nodes))
print("# Edges " + str(len(G.edges())) + "\n")

def print_corrected_read(cid, seq_to_ranges, nano_to_ranges, length):

    join = []
    for sid in seq_to_ranges:
        s, e, _ = seq_to_ranges[sid]
        join.append( (s,e) )
    for sid in nano_to_ranges:
        s, e = nano_to_ranges[sid]
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
        output.write( ">" + cid + "_" + str(i) + "\n")
        write_sequence(output, find_sequence_r(cid, max(cs, 200),min(ce, length - 200)))

## take arbitrary node of G (always the min for reproducability) as starting point of traversal with bredth-first search
## add passed nodes to subset until threshold size is reached
## the nodes of this connected subgraph are mapped all versus all
## all nodes from this subset that have all their neighbors in the set are scrubbed and then removed from G
## repeat until G is empty

# set size threshold of connected subset of nodes
bfs_subset = set() 
counter = 0

# start printing progress
print ('0/'+str(total_nodes))

while (len(G) != 0):
    if(len(bfs_subset) > 0):
        possible_starts = set(G.nodes()).difference(bfs_subset) # avoid to choose starting vertex that is already contained in subset
        
    else: possible_starts = set(G.nodes())
    start_node = min(possible_starts)
    
    #breadth-first-search
    bfs_edges = nx.bfs_edges(G, start_node, depth_limit = subset_size) 
    bfs_nodes = [start_node] + [v for u, v in bfs_edges]
        
    # add nodes starting from a random node to subset following bread-first-search until threshold is reached
    for node in bfs_nodes: 
        if(len(bfs_subset) >= subset_size):
            break
        else:
            bfs_subset.add(node)
    
    if(len(bfs_subset) < subset_size and len(G) > len(bfs_subset)): # the conncected component is smaller than threshold size -> so merge it to next one
        continue
    
    #copy the set and delete all nodes with neighbors outside the set
    bfs_subset_center = bfs_subset.copy()
    for u, v in nx.edge_boundary(G, bfs_subset):
        bfs_subset_center.discard(u)
    
    # map nodes of this subset against each other    
    # write sequences to file
    f1 = open(tmp_dir+"/temp_sequences.fa", "w")  
    for node in bfs_subset:
        f1.write(">"+node+"\n")
        write_sequence(f1, find_sequence(node))
    f1.close()

    os.system('minimap2 -x ava-ont ' + tmp_dir + '/temp_sequences.fa ' + tmp_dir + '/temp_sequences.fa > ' + tmp_dir + '/temp_pwa.paf 2> /dev/null')

    # read information of paf file and write to node dictionary seq_to_ranges
    for line in open(tmp_dir+'/temp_pwa.paf', "r"):

        linetemp = line.rstrip().split("\t")
        if len(linetemp) == 1:
            continue

        id_1, id_2 = linetemp[0], linetemp[5]  # remember about 0-based counting
        if id_1 == id_2:
            continue
  
        s_1, s_2 = int(linetemp[2]), int(linetemp[7])
        e_1, e_2 = int(linetemp[3]), int(linetemp[8])

        if e_1 - s_1 < 500:
            continue

        direction = linetemp[4]
        
        # write information for id_1 to dict 
        if id_2 not in G.nodes[id_1]['seq_to_ranges']:
            G.nodes[id_1]['seq_to_ranges'][id_2] = (s_1, e_1, direction)
        else:
            str_1 = G.nodes[id_1]['seq_to_ranges'][id_2][0]
            str_2 = G.nodes[id_1]['seq_to_ranges'][id_2][1]
            d = G.nodes[id_1]['seq_to_ranges'][id_2][2]
            if direction == d and ( abs( str_1 - e_1 ) < 500 or abs(s_1 - str_2) < 500):
                G.nodes[id_1]['seq_to_ranges'][id_2] = (min(s_1, str_1), max(e_1, str_2), direction)

        # write information for id_2 to dict
        if id_1 not in G.nodes[id_2]['seq_to_ranges']:
            G.nodes[id_2]['seq_to_ranges'][id_1] = (s_2, e_2, direction)
        else:
            str_1 = G.nodes[id_2]['seq_to_ranges'][id_1][0]
            str_2 = G.nodes[id_2]['seq_to_ranges'][id_1][1]
            d = G.nodes[id_2]['seq_to_ranges'][id_1][2]
            if direction == d and ( abs( str_1 - e_2 ) < 500 or abs(s_2 - str_2) < 500):
                G.nodes[id_2]['seq_to_ranges'][id_1] = (min(s_2, str_1), max(e_2, str_2), direction)

    # iterate over the nodes of the subset and write scrubbed sequences to output file
    for node in bfs_subset_center:
        print_corrected_read(node, G.nodes[node]['seq_to_ranges'], G.nodes[node]['illu_to_ranges'], G.nodes[node]['length'])
            
    # remove nodes from G
    G.remove_nodes_from(bfs_subset_center)
    subset_length = len(bfs_subset_center)
    bfs_subset.clear() 

    # print progress
    counter += subset_length
    print(str(counter)+'/'+str(total_nodes))
    
print("--- %s seconds ---" % (time.time() - start_time))
    
