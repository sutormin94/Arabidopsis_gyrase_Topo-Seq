###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

#Prepare motif for E. coli gyrase in Cfx conditions.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
from Bio import SeqIO

# Path to E. coli gyrase GSCs sequences.
GSCs_seq_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\\Results\OUTPUT_Motifs_visualization_sequences_extraction.py\\Cfx_sequences_under_GCSs_full.fasta"

# Path to output file with gyrase motif.
GCSs_motif_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\Motif_0.01\\E_coli_gyrase_Cfx_GC_pfm.txt"


#######
# Read file with GCSs sequences.
#######

def GCSs_seq(GSCs_seq_path):

    #Reads FASTA file with the reference genome. 
    seq_list=[]

    for record in SeqIO.parse(GSCs_seq_path, "fasta"):
        GCSs_id=record.id
        GSCs_seq=str(record.seq)

        seq_list.append(GSCs_seq)

    return seq_list


#######
#PFM construction.
#Scans sequences stack by columns, counts the number of particular letters.
#Returns a range of PFMs - "positional frequencies matrixes" .
#######

def make_PFM(seqs_list):

    matrix=[]
    template=seqs_list[0]
    for i in range(len(template)):
        column=[0, 0, 0, 0]
        for j in range(len(seqs_list)):
            if seqs_list[j][i] == str('A'):
                column[0] = column[0] + 1
            elif seqs_list[j][i] == str('T'):
                column[1] = column[1] + 1
            elif seqs_list[j][i] == str('G'):
                column[2] = column[2] + 1
            elif seqs_list[j][i] == str('C'):
                column[3] = column[3] + 1
        matrix.append(column)
    #Returns a range of PFMs.
    GC_percent = []
    GT_percent = []
    CT_percent = []
    A_percent = []
    T_percent = []
    G_percent = []
    C_percent = []
    for i in range(len(matrix)):
        GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
                    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GC_percent.append(GC)
        GT_percent.append(GT)
        CT_percent.append(CT)
        A_percent.append(A)
        T_percent.append(T)
        G_percent.append(G)
        C_percent.append(C)

    return {'Num_seqs': len(seqs_list), 'A': A_percent, 'T': T_percent, 'G': G_percent, 'C': C_percent, 'CT': CT_percent, 'GT': GT_percent, 'GC': GC_percent}


#######
#Writes PFM data to file.
#######

def write_motif(ar, filepath, coord_shift):

    fileout=open(filepath, 'w')
    fileout.write("#X\tY\n")
    for i in range(len(ar)):
        fileout.write(str((-coord_shift/2)+1+i) + '\t' + str(ar[i])+'\n')
    fileout.close()

    return



def wrapper_func(sequences_list_path, output_path):
    
    win_width=170
    PFM_type='GC'
    
    sequences_list=GCSs_seq(sequences_list_path)
    PFMs=make_PFM(sequences_list)
    write_motif(PFMs[PFM_type], output_path, win_width)
    
    return

wrapper_func(GSCs_seq_path, GCSs_motif_path)