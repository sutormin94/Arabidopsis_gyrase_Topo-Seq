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
PWD_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\"

# Path to fasta with all chromosomes.
Fasta_path=os.path.join(PWD_path, "Reference_genome", "GCA_000001735.2_TAIR10.1_genomic_NUMT_masked.fna")

# Path to gff with all chromosomes.
GFF_path=os.path.join(PWD_path, "Reference_genome", "GCA_000001735.2_TAIR10.1_genomic.gff")

# Path to WIG files with all chromosomes.
WIG_path_dict={"Topo_Seq_1_S5_edt_for_rev_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_for_rev_depth.wig"),
               "Topo_Seq_1_S5_edt_forward_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_forward_depth.wig"),
               "Topo_Seq_1_S5_edt_reverse_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_reverse_depth.wig"),
               "Topo_Seq_1_S5_edt_N3E_F"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_N3E_F.wig"),
               "Topo_Seq_1_S5_edt_N3E_R"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_N3E_R.wig"),
               "Topo_Seq_1_S5_edt_N5E_F"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_N5E_F.wig"),
               "Topo_Seq_1_S5_edt_N5E_R"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_1_S5_edt_N5E_R.wig"),
               "Topo_Seq_2_S6_edt_for_rev_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_for_rev_depth.wig"),
               "Topo_Seq_2_S6_edt_forward_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_forward_depth.wig"),
               "Topo_Seq_2_S6_edt_reverse_depth"      : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_reverse_depth.wig"),
               "Topo_Seq_2_S6_edt_N3E_F"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_N3E_F.wig"),
               "Topo_Seq_2_S6_edt_N3E_R"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_N3E_R.wig"),
               "Topo_Seq_2_S6_edt_N5E_F"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_N5E_F.wig"),
               "Topo_Seq_2_S6_edt_N5E_R"              : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_2_S6_edt_N5E_R.wig"),
               "Topo_Seq_Mock_1_S7_edt_for_rev_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_for_rev_depth.wig"),
               "Topo_Seq_Mock_1_S7_edt_forward_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_forward_depth.wig"),
               "Topo_Seq_Mock_1_S7_edt_reverse_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_reverse_depth.wig"),
               "Topo_Seq_Mock_1_S7_edt_N3E_F"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_N3E_F.wig"),
               "Topo_Seq_Mock_1_S7_edt_N3E_R"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_N3E_R.wig"),
               "Topo_Seq_Mock_1_S7_edt_N5E_F"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_N5E_F.wig"),
               "Topo_Seq_Mock_1_S7_edt_N5E_R"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_1_S7_edt_N5E_R.wig"),
               "Topo_Seq_Mock_2_S8_edt_for_rev_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_for_rev_depth.wig"),
               "Topo_Seq_Mock_2_S8_edt_forward_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_forward_depth.wig"),
               "Topo_Seq_Mock_2_S8_edt_reverse_depth" : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_reverse_depth.wig"),
               "Topo_Seq_Mock_2_S8_edt_N3E_F"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_N3E_F.wig"),
               "Topo_Seq_Mock_2_S8_edt_N3E_R"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_N3E_R.wig"),
               "Topo_Seq_Mock_2_S8_edt_N5E_F"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_N5E_F.wig"),
               "Topo_Seq_Mock_2_S8_edt_N5E_R"         : os.path.join(PWD_path, "NUMT_masked", "WIG", "Topo_Seq_Mock_2_S8_edt_N5E_R.wig"),
               "ChIP_Seq_2_S2"                        : os.path.join(PWD_path, "NUMT_masked", "Gyrase_ChIP_Seq", "WIG", "TAIR10.1_ChIP_Seq_NUMT_masked_2_S2.wig"),
               "ChIP_Seq_Mock_2_S4"                   : os.path.join(PWD_path, "NUMT_masked", "Gyrase_ChIP_Seq", "WIG", "TAIR10.1_ChIP_Seq_NUMT_masked_Mock_2_S4.wig"),
               "ChIP_Seq_FE_2"                        : os.path.join(PWD_path, "NUMT_masked", "Gyrase_ChIP_Seq", "WIG", "TAIR10.1_ChIP_Seq_NUMT_masked_FE_2.wig"),
               }

# Path to BroadPeak files with GSCs coordinates.
GCSs_path_dict={"Trusted_GCSs_0.01" : os.path.join(PWD_path, "NUMT_masked", "Trusted_GSCs", "Trusted_GCSs_0.01", "Trusted_GCSs_0.01_all_chromosomes.BroadPeak"),
                "Trusted_GCSs_0.05" : os.path.join(PWD_path, "NUMT_masked", "Trusted_GSCs", "Trusted_GCSs_0.05", "Trusted_GCSs_0.05_all_chromosomes.BroadPeak"),
                }

# List of chromosomes to be extracted.
Chrom_ar=["AP000423.1", "BK010421.1"]

# Path for output.
Outpath=os.path.join(PWD_path, "NUMT_masked", "Separate_chromosomes_data")
if os.path.isdir(Outpath)==False:
    os.mkdir(Outpath)



def read_split_fasta(fasta_path, chrom_ar, outpath):
    
    #Reads FASTA file with the reference genome. 
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        chr_id=record.id
        chr_seq=str(record.seq)
        
        if chr_id in chrom_ar:
            if os.path.isdir(os.path.join(outpath, chr_id))==False:
                os.mkdir(os.path.join(outpath, chr_id))
                
            print(f'Extracting chromosome {chr_id} sequence, length {len(chr_seq)} nt')  
                
            fileout=open(os.path.join(outpath, chr_id, f'{chr_id}.fasta'), 'w')
            fileout.write(f'>{chr_id}\n{chr_seq}')
            fileout.close()

    return


def read_split_gff(gff_path, chrom_ar, outpath):
    
    for extract_chrom in chrom_ar:
        
        if os.path.isdir(os.path.join(outpath, extract_chrom))==False:
            os.mkdir(os.path.join(outpath, extract_chrom)) 
            
        print(f'Extracting chromosome {extract_chrom} gff annotation') 
    
        filein=open(gff_path, 'r')
        fileout=open(os.path.join(outpath, extract_chrom, f'{extract_chrom}.gff'), 'w')
        
        for line in filein:
            if extract_chrom in line:
                fileout.write(line)
        
        filein.close()
        fileout.close()
    
    return


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


def read_split_BroadPeak(gcss_path_dict, chrom_ar, outpath):
    
    for broadpeak_name, broadpeak_path in gcss_path_dict.items():
        
        print(f'Processing {broadpeak_name} BroadPeak annotation')
    
        for extract_chrom in chrom_ar:
            
            if os.path.isdir(os.path.join(outpath, extract_chrom))==False:
                os.mkdir(os.path.join(outpath, extract_chrom)) 
                
            print(f'Extracting chromosome {extract_chrom} BroadPeak annotation') 
        
            filein=open(broadpeak_path, 'r')
            fileout=open(os.path.join(outpath, extract_chrom, f'{extract_chrom}_{broadpeak_name}.BroadPeak'), 'w')
            
            for line in filein:
                if extract_chrom in line:
                    fileout.write(line)
            
            filein.close()
            fileout.close()    
    
    return


def wrapper_function(fasta_path, gff_path, wig_path_dict, gcss_path_dict, chrom_ar, outpath):
    
    # Read and split fasta file.
    read_split_fasta(fasta_path, chrom_ar, outpath)
    
    # Read and split gff file.
    read_split_gff(gff_path, chrom_ar, outpath)
    
    # Read and split wig files.
    read_split_wig(wig_path_dict, chrom_ar, outpath)
    
    # Read and split BroadPeak files.
    read_split_BroadPeak(gcss_path_dict, chrom_ar, outpath)
    
    
    return

wrapper_function(Fasta_path, GFF_path, WIG_path_dict, GCSs_path_dict, Chrom_ar, Outpath)