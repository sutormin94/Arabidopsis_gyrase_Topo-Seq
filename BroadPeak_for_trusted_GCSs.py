###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis Topo-Seq analysis##

# Prepares BroadPeak files containing trusted GCSs coordinates for multiple chromosomes
# from lists of trusted GCSs made for individual chromosomes.
###############################################

#######
#Packages to be imported.
#######

import os

# Path to working directory.
Output_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\\Trusted_GCSs_0.05\\"

# Dictionary of files with trusted GCSs for individual chromosomes.
GCSs_path_dict={"CP002684.1" : os.path.join(Output_path, 'CP002684.1_Chr1_GSCs_trusted_GCSs.txt'),
                "CP002685.1" : os.path.join(Output_path, 'CP002685.1_Chr2_GSCs_trusted_GCSs.txt'),
                "CP002686.1" : os.path.join(Output_path, 'CP002686.1_Chr3_GSCs_trusted_GCSs.txt'),
                "CP002687.1" : os.path.join(Output_path, 'CP002687.1_Chr4_GSCs_trusted_GCSs.txt'),
                "CP002688.1" : os.path.join(Output_path, 'CP002688.1_Chr5_GSCs_trusted_GCSs.txt'),
                "AP000423.1" : os.path.join(Output_path, 'AP000423.1_Chloroplast_GSCs_trusted_GCSs.txt'),
                "BK010421.1" : os.path.join(Output_path, 'BK010421.1_Mitochondria_GSCs_trusted_GCSs.txt'),
                }

# Dataset name.
Dataset_name="Trusted_GCSs_0.05_all_chromosomes.BroadPeak"


def wrapper_function(dataset_name, gcss_path_dict, output_path):
    
    output_file_path=os.path.join(output_path, dataset_name)
    fileout=open(output_file_path, 'w')
    
    GSCs_ID=0
    
    for chr_name, file_path in gcss_path_dict.items():
        
        filein=open(file_path, 'r')
        
        for line in filein:
            
            line=line.rstrip().split('\t')
            
            if line[0] not in ['GCSs_coordinate']:
                
                coord_start=int(line[0])
                coord_end=coord_start+4
                N3E=float(line[1])
                
                fileout.write(f'{chr_name}\t{coord_start}\t{coord_end}\tGCS_{GSCs_ID}\t{N3E}\t.\t-1\t-1\t-1\n')
                
                GSCs_ID+=1
                
        filein.close()
        
    fileout.close()
    
    return

wrapper_function(Dataset_name, GCSs_path_dict, Output_path)