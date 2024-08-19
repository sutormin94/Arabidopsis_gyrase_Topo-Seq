###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis Topo-Seq analysis##

# Extract specified chromosomes data from fasta/gff/wig/BroadPeak files.
###############################################

#######
#Packages to be imported.
#######

import os
from Bio import SeqIO


# Path to working directory.
PWD_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\Mito_RNA_Seq_Garcia_2021\WIG\\"

# Path to WIG files with all chromosomes.
WIG_path_dict={"MtRNA_Seq_1_ERR1665215"   : os.path.join(PWD_path, "MtRNA_Seq_1_ERR1665215.wig"),
               "MtRNA_Seq_2_ERR1665219"   : os.path.join(PWD_path, "MtRNA_Seq_2_ERR1665219.wig"),
               "MtRNA_Seq_3_ERR1665220"   : os.path.join(PWD_path, "MtRNA_Seq_3_ERR1665220.wig"),
               }

# List of chromosomes to be extracted.
Chrom_ar=["BK010421.1"]

# Path for output.
Outpath=os.path.join(PWD_path)
if os.path.isdir(Outpath)==False:
    os.mkdir(Outpath)


#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    
    wigin=open(wigfile, 'r')
    NE_values_dict={}
    NE_values=[]
    
    for line in wigin:
        line=line.rstrip().split(' ')
        
        if line[0]!='track':
        
            if line[0]=='fixedStep':
                
                if len(NE_values)>0:
                    
                    NE_values_dict[chr_id]=NE_values
                
                chr_id=line[1].split('=')[1]
                NE_values=[]               
                
            else:
                
                NE_values.append(float(line[0]))
            
    NE_values_dict[chr_id]=NE_values
    
    wigin.close()
    
    return NE_values_dict


def read_split_wig(wig_path_dict, chrom_ar, outpath):
    
    for wig_track_name, wig_path in wig_path_dict.items():
        
        print(f'Processing {wig_track_name} wig tracks')
        
        Chrom_tracks_dict=wig_parsing(wig_path)
        
        for extract_chrom in chrom_ar:
            
            if extract_chrom in Chrom_tracks_dict:
                
                print(f'Extracting chromosome {extract_chrom} wig track')
                
                if os.path.isdir(os.path.join(outpath, extract_chrom))==False:
                    os.mkdir(os.path.join(outpath, extract_chrom))
                    
                fileout=open(os.path.join(outpath, extract_chrom, f'{extract_chrom}_{wig_track_name}.wig'), 'w')
                
                fileout.write(f'track type=wiggle_0 name={wig_track_name} autoScale=off viewLimits=0.0:25.0\nfixedStep chrom={extract_chrom} start=1 step=1\n')
                
                for value in Chrom_tracks_dict[extract_chrom]:
                    
                    fileout.write(f'{value}\n')
                
                fileout.close()
                
            else:
                
                print(f'Chromosome {extract_chrom} not found in wig track {wig_path}')
    
    return





def wrapper_function(wig_path_dict, chrom_ar, outpath):
    
    # Read and split wig files.
    read_split_wig(wig_path_dict, chrom_ar, outpath)
    
    return

wrapper_function(WIG_path_dict, Chrom_ar, Outpath)