###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase analysis##

####
#Converts wig file with some genomic feature to Circos-compatible format (bed-like).
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import scipy
import os

# Path to working directory.
PWD_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\Separate_chromosomes_data\\"

# Input files.
Input_path=os.path.join(PWD_path, 'BK010421.1', 'BK010421.1_Topo_Seq_1_S5_edt_reverse_depth.wig')

# Bin width, nt.
Bin_width=2

# Scale factor.
Scale_factor=-1

# Chromosome ID.
Chromosome_id='BK010421.1'

# Output file.
Output_path=os.path.join(PWD_path, 'BK010421.1', 'BK010421.1_Topo_Seq_1_S5_edt_reverse_depth_bin_2_negative.txt')


#######
##Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(f'Track length: {len(NE_values)} bp')
    return NE_values

#######
##Read wig, bin, convert to Circus format.
#######

def read_bin_write(inpath, bwidth, scale_factor, chrid, outpath):
    #Read wig data.
    wig_data=wig_parsing(inpath)
    
    #Bin wig data.
    wig_binned=[]
    start=0
    while start+bwidth<len(wig_data):
        wig_binned.append(scale_factor*np.mean(wig_data[start:(start+bwidth)]))
        start+=bwidth
    wig_binned.append(scale_factor*np.mean(wig_data[start:]))
    
    
    #Write binned data in Circos format.
    outfile=open(outpath, 'w')
    for i in range(len(wig_binned)):
        if ((i+1)*bwidth)<len(wig_data):
            outfile.write(f'{chrid} {i*bwidth} {(i+1)*bwidth} {wig_binned[i]}\n')
        else:
            outfile.write(f'{chrid} {i*bwidth} {len(wig_data)} {wig_binned[i]}\n')
    
    outfile.close()
    return

read_bin_write(Input_path, Bin_width, Scale_factor, Chromosome_id, Output_path)