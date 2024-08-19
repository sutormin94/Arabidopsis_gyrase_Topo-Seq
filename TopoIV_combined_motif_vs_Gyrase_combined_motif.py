###############################################
##Dmitry Sutormin, 2024##
##Topo-Seq analysis##

#The script takes several PFM (e.g. Gyrase from E. coli and A. thaliana) and plots them together. 
###############################################

#######
#Packages to be imported.
#######

import os
import scipy
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet  import IUPAC
import matplotlib.pyplot as plt

#######
#Variables to be defined.
#######

print('Variables to be defined:')

PWD="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\\"

#Input: Motifs, two-column TAB.
path_to_combined_motifs={'E. coli' :            PWD + "Motif_0.05\E_coli_gyrase_Cfx_GC_pfm.txt", 
                         'A. thaliana Chloro' : PWD + "Motif_0.05\AP000423.1_GC_pfm.txt",
                         'A. thaliana Mito' :   PWD + "Motif_0.05\BK010421.1_GC_pfm.txt",
                         }

#Output: prefix of the output path.
Output_data_prefix=PWD + "Motif_0.05\\"
if not os.path.exists(Output_data_prefix):
    os.makedirs(Output_data_prefix)
    
    
###############################################
#Reading the motif files and visualization.
###############################################

#######
#Motif parsing.
#######

def read_motif(inpath):
    filein=open(inpath, 'r')
    Coords=[]
    GC=[]
    for line in filein:
        if line[0]!='#':
            line=line.rstrip().split('\t')
            Coords.append(int(float(line[0])))
            GC.append(float(line[1]))
    filein.close()
    return Coords, GC


#######
#Plots combined motifs.
#######

def Plotting(Motifs_dict, Coordinates, matrix_type, outpath):
    
    yticks=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    xticks=[-63, -11, 0, 14, 66]
    xticks_minor=[-80, -70, -60, -50, -40, -30, -20, -10, 10, 20, 30, 40, 50, 60, 70, 80]
    x_axis=Coordinates
    ax_range=[min(x_axis), max(x_axis), 0.2, 0.9]
    Colors_palette=[['#ff9fda', '#bf76a1'], ['#ffb486', '#cc8e67'], ['#8bdaff', '#67a0ba']]
    plt.figure(figsize=(16, 6), dpi=100)
    plot1=plt.subplot()
    i=0
    
    for motif_name, motif in Motifs_dict.items():
        
        plot1.plot(x_axis, motif, color=Colors_palette[i][0], linewidth=7, alpha=0.6, zorder=1)
        plot1.plot(x_axis, motif, color=Colors_palette[i][1], linewidth=3, alpha=0.6, zorder=3, label=motif_name)
        plot1.plot(x_axis, motif, 'o', fillstyle='none', color=Colors_palette[i][0], markeredgecolor=Colors_palette[i][1], markersize=7, alpha=0.6, zorder=2)               
        plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
        plot1.axis(ax_range)
        plot1.set_xlim(min(x_axis), max(x_axis))
        
        plot1.set_xticks(xticks, minor=False)
        plot1.tick_params(axis='both', which='major', labelsize=35)
        plot1.xaxis.grid(True, which='major', linewidth=0.5, linestyle='--', alpha=1)     
        
        plot1.set_xticks(xticks_minor, minor=True)
        plot1.tick_params(axis='x', which='minor', labelsize=20)
        
        plot1.set_yticks(yticks, minor=False)
        
        #plot1.set_xticks([0], minor=True)
        #plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])), minor=False)
        
        plot1.set_xlabel('Position, nt', size=35)
        plot1.set_ylabel(str(matrix_type + '%'), size=35)
        i+=1
        
    plt.legend(fontsize=25, frameon=False, markerscale=3, handlelength=0.5, bbox_to_anchor=(1, 1)) # ((0.7, 0.6) and labelspacing=0.4) or ((1, 1) and no labelspacing=0.4)
    plt.tight_layout()
    plt.savefig(outpath + 'Combined_motif_Gyrase_A_thaliana_and_E_coli.png', dpi=400, figsize=(16, 6))
    plt.savefig(outpath + 'Combined_motif_Gyrase_A_thaliana_and_E_coli.svg', dpi=400, figsize=(16, 6))
    plt.close()
    
    return

#######
#Wrapper function.
#######

def wrapper_motif(path_motif_dict, matrix_type, outpath):
    
    Motifs_dict={}
    for name, inpath in path_motif_dict.items():
        Coordinates, Motif=read_motif(inpath)
        Motifs_dict[name]=Motif
    Plotting(Motifs_dict, Coordinates, matrix_type, outpath)
    
    return

wrapper_motif(path_to_combined_motifs, 'GC', Output_data_prefix)