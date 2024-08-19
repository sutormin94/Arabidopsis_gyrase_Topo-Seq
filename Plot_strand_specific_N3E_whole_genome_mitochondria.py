###############################################
##Dmitry Sutormin, 2024##
##Topo-Seq analysis. Plotting of strand-specific data.##

#The script takes strand-specific tracks (e.g., N3E) and plots a specified genomic region.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.patches as patches
import numpy as np


#######
#Variables to be defined.
#######

# Path to the working directory
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\"

# Dataset name.
Dataset_name="BK010421.1_replicate_1_cp_ori_masked"

# N3E data.
N3E_tracks={'IP_N3E_F'  : os.path.join(PWD, "WIG", "BK010421.1", "BK010421.1_Topo_Seq_1_S5_edt_N3E_F.wig"), 
            'IP_N3E_R'  : os.path.join(PWD, "WIG", "BK010421.1", "BK010421.1_Topo_Seq_1_S5_edt_N3E_F.wig"),             
           }

# Mock data.
Mock_track_path=os.path.join(PWD, "WIG", "BK010421.1", "BK010421.1_Topo_Seq_Mock_1_S7_edt_for_rev_depth_cp_ori_masked.wig")

# Path to file with GCSs coordinates.
GCS_coord_path=os.path.join(PWD, "Separate_chromosomes_data", "BK010421.1", "BK010421.1_Trusted_GCSs_0.01.BroadPeak")

# Genomic region to plot.
Region_to_plot=[0, 367808]




#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    Total_NE=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
            Total_NE+=float(line[0])
    print('Total number of ends: ' + str(Total_NE))
    wigin.close()
    return NE_values, Total_NE


#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar


#######
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
    
    if i<0: #coordinate is out of the left genome border (start)
        j=len(ends)+i
    elif i>=len(ends): #coordinate is out of the right genome border (end)
        j=i-len(ends)
    else: #coordinate is within the genome borders
        check_in_del=0
        for dl in deletions: #check if coordinate falls into deletion
            if dl[1]>=i>=dl[0]:
                j=dl[1]-dl[0]+i+1
                check_in_del=1
        if check_in_del==0:
            j=i
            
    return ends[j]


#######
#Returns smoothed N3/5E tracks.
#Smoothing using sliding window.
#######

def Smoothing(ends, deletions):
    
    smoothed=[]
    #Calculating the value for the first genome position
    mean=0.0
    window=10000
    window_float=float(window)
    for i in range(-window, window):
        mean=mean + get_value(i, ends, deletions)
    mean=mean/(2*window_float)
    smoothed.append(mean)
    #Calculating values for the part of the genome remains
    for i in range(1, len(ends)):
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
        
    return smoothed


######
#Plots signal over the genomic interval.
#######

def plot_signal(dataset_name, Wig_tracks_dict, Mock_track_smoothed, TCSs_coords, region_to_plot, path_out):
    
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    
    #Get left-most TCSs positions.
    TSCs_ar=[]
    for TCS in TCSs_coords:
        TCS_left=TCS[0]
        TSCs_ar.append(TCS_left)    
    
    #Identify max value.
    wig_max=0
    for wig_name, wig_data in Wig_tracks_dict.items():
        if max(wig_data)>wig_max:
            wig_max=max(wig_data)
    wig_max=wig_max/2
    
    #Plotting the distribution of the signal around the genome for IPed and mock samples.
    Ori_1=272000
    Ori_2=333500
    CpOri_1=43000
    CpOri_2=129900
    CpOri_3=165200
    Ori_coords=[Ori_1, Ori_2, CpOri_1, CpOri_2, CpOri_3]
    xcoord=np.arange(0,len(wig_data))
    fig, plots=plt.subplots(3, figsize=(8, 5), dpi=100)
    
    for wig_name, wig_data in Wig_tracks_dict.items():
        
        if ('IP_' in wig_name) and ('_F' in wig_name):
            plots[0].plot(xcoord, wig_data, '-', label=wig_name, color='black', linewidth=0.5)
            plots[0].set_ylabel('N3E F', size=17)
            plots[0].set_ylim([0,wig_max])
            plots[0].set_xlim(region_to_plot)
            plots[0].plot(xcoord, Mock_track_smoothed, '-', label=f'Mock, smoothed', color='orange', linewidth=3)
            plots[0].vlines(Ori_coords, 0, wig_max, color='red', linestyles="dashed", linewidth=1)
            plots[0].spines["top"].set_visible(False)
            plots[0].spines["right"].set_visible(False) 
            plots[0].legend(loc='upper right', frameon=False)
            
        elif ('IP_' in wig_name) and ('_R' in wig_name):
            plots[1].plot(xcoord, wig_data, '-', label=wig_name, color='black', linewidth=0.5)
            plots[1].vlines(Ori_coords, 0, wig_max, color='red', linestyles="dashed", linewidth=1)
            plots[1].set_ylabel('N3E R', size=17)
            plots[1].set_xlim(region_to_plot)
            plots[1].set_ylim([0,wig_max])
            plots[1].spines["top"].set_visible(False)
            plots[1].spines["right"].set_visible(False)   
            plots[1].legend(loc='upper right', frameon=False)
    
    plots[2].set_xlim(region_to_plot)
    plots[2].set_ylim([-0.03, 0.03])
    plots[2].scatter(TSCs_ar, [0]*len(TSCs_ar), s=2, marker='o', c='black', alpha=0.1, label='TCSs')
    plots[2].vlines(Ori_coords, -0.03, 0.03, color='red', linestyles="dashed", linewidth=1)
    
    plt.show()
    plt.savefig(os.path.join(path_out, f'{dataset_name}_strand_specific_signal.png'), dpi=300, figsize=(8, 5))
    plt.savefig(os.path.join(path_out, f'{dataset_name}_strand_specific_signal.svg'), dpi=300, figsize=(8, 5))
    plt.close()
    
    return


def wrapper_func(dataset_name, pwd, s_s_tracks, mock_track_path, gcs_coord_path, region_to_plot):
    
    # Read N3E wig tracks.
    Wig_tracks_dict={}
    for wig_name, wig_path in s_s_tracks.items():
        Wig_tracks_dict[wig_name]=wig_parsing(wig_path)[0]
        
    # Read Mock wig track.
    Deletions=[]
    Mock_track=wig_parsing(mock_track_path)[0]
    Mock_track_smoothed=Smoothing(Mock_track, Deletions)
    
    #Parsing TCSs coords.
    TCSs_coords=deletions_info(gcs_coord_path)    
    
    # Plot wig tracks.
    plot_signal(dataset_name, Wig_tracks_dict, Mock_track_smoothed, TCSs_coords, region_to_plot, pwd)
    
    return


wrapper_func(Dataset_name, PWD, N3E_tracks, Mock_track_path, GCS_coord_path, Region_to_plot)