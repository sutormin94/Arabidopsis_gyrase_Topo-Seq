###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over transcription units (TUs). Plots this information.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np


# Path to working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\WIG\BK010421.1\\"

# Track to mask.
Wig_track_path=os.path.join(PWD, "BK010421.1_Topo_Seq_Mock_2_S8_edt_for_rev_depth.wig")

# Regions to mask.
Regions_to_mask=[[42800, 43400], [129600, 130100], [165000, 165250]]

# Output path.
Wig_track_outpath=os.path.join(PWD, "BK010421.1_Topo_Seq_Mock_2_S8_edt_for_rev_depth_cp_ori_masked.wig")




#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    
    header_ar=[]
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
        else:
            header_ar.append(line)
            
    wigin.close()
    
    return header_ar, NE_values


def mask_wig(NE_values, regions_to_mask):
    
    wig_track_masked=[]
    
    mean_signal=np.mean(NE_values)
    
    for mask_reg in regions_to_mask:
        
        reg_start=mask_reg[0]
        reg_end=mask_reg[1]
        
        NE_values[reg_start:reg_end]=[mean_signal]*(reg_end-reg_start)
    
    return NE_values


def write_wig(header_ar, wig_track_masked, wig_track_outpath):
    
    fileout=open(wig_track_outpath, 'w')
    
    for header_line_ar in header_ar:
        
        for word in header_line_ar:
            
            fileout.write(f'{word} ')
        fileout.write(f'\n')
        
    for value in wig_track_masked:
        
        fileout.write(f'{value}\n')
    
    fileout.close()
    
    return


def wrapper_func(wig_track_path, regions_to_mask, wig_track_outpath):
    
    # Read wig.
    header_ar, wig_track=wig_parsing(wig_track_path)
    
    # Mask regions.
    wig_track_masked=mask_wig(wig_track, regions_to_mask)
    
    # Write masked wig.
    write_wig(header_ar, wig_track_masked, wig_track_outpath)
    
    return

wrapper_func(Wig_track_path, Regions_to_mask, Wig_track_outpath)