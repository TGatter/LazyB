#!/usr/bin/env python

import numpy as np
import sys
from Bio import SeqIO

paf_file = sys.argv[1]
unitigs = sys.argv[2]
output_stats = sys.argv[3]
output_file = sys.argv[4]

unitigs_corrected = open(output_file, 'w')
stats_file = open(output_stats, 'a')

# open up sequence map for illumina
illumina_dict = SeqIO.index_db(unitigs + ".idx", unitigs, "fasta")


def find_sequence_illumina(seqid, start, end):
    return illumina_dict[seqid].seq[start:end + 1]


def write_sequence(f, seq):
    for i in range(0, len(seq), 60):
        f.write(str(seq[i:i + 60]) + "\n")


############ read input paf file #######################################################################################
mapping_freq = {}
nano_lengths = {}

with open(paf_file, 'r') as input_file:
    current_id = ''
    current_length = 0
    current_nano_ids = []
    current_profile = []
    block_counter = 0

    for line in input_file:
        split_line = line.rstrip().split('\t')
        illu_id = split_line[0]
        illu_length = int(split_line[1])
        illuminaStart = int(split_line[2])
        illuminaEnd = int(split_line[3]) - 1
        nano_id = split_line[5]
        nano_length = split_line[6]

        if nano_id not in nano_lengths:
            nano_lengths[nano_id] = int(nano_length)

        if current_id == illu_id:  # still in the same block
            if nano_id not in current_nano_ids:
                block_counter += 1
                current_nano_ids.append(nano_id)
                for pos in range(illuminaStart, illuminaEnd + 1):
                    current_profile[pos] += 1

        else:  # new block
            if current_id != '':  # close previous block, avoid beginning
                mapping_freq[current_id] = [block_counter, current_length, max(current_profile)]

            # start new block
            current_profile = [0] * illu_length
            for pos in range(illuminaStart, illuminaEnd + 1):
                current_profile[pos] += 1

            block_counter = 1
            current_id = illu_id
            current_length = illu_length
            del current_nano_ids[:]
            current_nano_ids.append(nano_id)

    mapping_freq[current_id] = [block_counter, current_length, max(current_profile)]  # last block
#print("finished reading input file")


############# define filter cut-off ####################################################################################

def filter_cutoff(mapping_dict):
    blocks_maxcov = [x[2] for x in list(mapping_dict.values())]
    # define cut-off 
    q1 = np.percentile(blocks_maxcov, 25)
    q3 = np.percentile(blocks_maxcov, 75)
    iqr = q3 - q1
    upper_whisker = q3 + 1.5 * iqr
    stats_file.write(">>> unitig filter \n")
    stats_file.write("upper_outlier: " + str(upper_whisker)+'\n')
    stats_file.write("Q3: " + str(q3)+'\n')

    return upper_whisker, blocks_maxcov, q3


cutoff, blocks_maxcov, q3 = filter_cutoff(mapping_freq)


########### evaluate and cut outliers  #################################################################################

def cut_peaks(unitig_id, coverage_profile, file):
    subrange = (0, 0)
    below_cutoff = False
    unitig_fragments = []
    for i, cov in enumerate(coverage_profile):
        if cov <= q3:                 # set q3 as cutoff
            if not below_cutoff:      # start new subrange below cutoff
                subrange = (i, i)
            else:                     # prolonging subrange
                subrange = (subrange[0], i)
            below_cutoff = True

        else:
            if below_cutoff and subrange[1] - subrange[0] + 1 >= 500:
                file.write('>' + unitig_id + '_' + str(len(unitig_fragments)) + ' ' + str(subrange[1] - subrange[0] + 1)
                           + ' ' + str(subrange[0]) + ' ' + str(subrange[1]) + '\n')
                write_sequence(file, find_sequence_illumina(unitig_id, subrange[0], subrange[1]))
                unitig_fragments.append(subrange)
                subrange = (0, 0)
            below_cutoff = False
    if below_cutoff and subrange[1] - subrange[0] + 1 >= 500:
        file.write('>' + unitig_id + '_' + str(len(unitig_fragments)) + ' ' + str(subrange[1] - subrange[0] + 1)
                   + ' ' + str(subrange[0]) + ' ' + str(subrange[1]) + '\n')
        write_sequence(file, find_sequence_illumina(unitig_id, subrange[0], subrange[1]))
        unitig_fragments.append(subrange)
    return unitig_fragments


########### detect outlier and write new unitig file   ###########################################################

with open(paf_file, 'r') as input_file:
    current_id = ''
    current_profile = []
    all_counter = 0
    outlier_counter = 0
    rescued_outlier = 0
    total_rescued_length=0
    for line in input_file:
        split_line = line.rstrip().split('\t')
        illu_id = split_line[0]
        illu_length = int(split_line[1])
        illuminaStart = int(split_line[2])
        illuminaEnd = int(split_line[3]) - 1

        if current_id == illu_id:  # still in the same block
            if mapping_freq[illu_id][2] > cutoff:  # outlier unitig block
                for pos in range(illuminaStart, illuminaEnd + 1):
                    current_profile[pos] += 1
            else:  # normal unitig block
                continue

        else:  # new block
            all_counter += 1
            if current_id != '' and mapping_freq[current_id][2] > cutoff:  # close previous block if outlier
                fragments = cut_peaks(current_id, current_profile, unitigs_corrected)
                if len(fragments) > 0:
                    rescued_outlier += 1
                    for entry in fragments:
                        total_rescued_length+=(entry[1] - entry[0] + 1)
                outlier_counter += 1

            if mapping_freq[illu_id][2] > cutoff:  # new outlier block
                current_profile = [0] * illu_length
                for pos in range(illuminaStart, illuminaEnd + 1):
                    current_profile[pos] += 1
                current_id = illu_id

            else:  # new normal unitig
                unitigs_corrected.write('>' + illumina_dict[illu_id].description + '\n')
                write_sequence(unitigs_corrected, illumina_dict[illu_id].seq)
                current_id = illu_id

    if current_id != '' and mapping_freq[current_id][2] > cutoff:  # close last block if outlier
        fragments = cut_peaks(current_id, current_profile, unitigs_corrected)
        if len(fragments) > 0:
            rescued_outlier += 1
        outlier_counter += 1

stats_file.write("#all unitigs: " + str(all_counter)+'\n')
stats_file.write("#outliers: " + str(outlier_counter)+'\n')
stats_file.write("#rescued outliers: " + str(rescued_outlier)+'\n')

unitigs_corrected.close()
stats_file.close()
