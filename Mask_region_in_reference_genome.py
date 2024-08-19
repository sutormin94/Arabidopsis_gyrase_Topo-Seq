###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

# Replaces a specified region of a specified chromosome with Ns.
# Extracts resultant sequence of the edited chromosome in a separate file.
# Returns a full reference genome (multiple chrommosomes) with edited chromosome specified.
###############################################

import Bio
from Bio import SeqIO

# ID of a chromosome to edit.
Chromosome_name="CP002685.1"
# Region start site to mask with Ns.
Reg_start=3239165
# Region end site to mask with Ns.
Reg_end=3509770

# Path to initial ref genome.
Old_ref_genome_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\Reference_genome\\GCA_000001735.2_TAIR10.1_genomic.fna"
# Output path to resultant edited genome.
New_ref_genome_path="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\Reference_genome\Edited_genome\\GCA_000001735.2_TAIR10.1_genomic_NUMT_masked.fna"
# Output path to fasta with a sequence of edited chromosome.
Edit_chrom_path=f"C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\Reference_genome\Edited_genome\{Chromosome_name}_NUMT_masked.fna"

def mask_NUMT(Old_ref_genome_path, New_ref_genome_path, Chromosome_name, Reg_start, Reg_end, Edit_chrom_path):
    
    genome_in=open(Old_ref_genome_path, 'r')
    genome_out=open(New_ref_genome_path, 'w')
    chrom_out=open(Edit_chrom_path, 'w')
    
    for record in SeqIO.parse(genome_in, "fasta"):
        sequence_name=record.description
        sequence=str(record.seq)
        init_seq_len=len(sequence)
        
        if Chromosome_name in sequence_name:
        
            print(sequence_name)
            
            Subst_seq="N"*(Reg_end-Reg_start)
            
            sequence=sequence[:Reg_start] + Subst_seq + sequence[Reg_end:]
            
            chrom_out.write(f'>{sequence_name}\n{sequence}\n')
            
        fin_seq_len=len(sequence)
        print(f'Init seq len: {init_seq_len}; final seq len: {fin_seq_len}')
            
        genome_out.write(f'>{sequence_name}\n{sequence}\n')
    
    
    genome_in.close()
    genome_out.close()
    chrom_out.close()
    
    return

mask_NUMT(Old_ref_genome_path, New_ref_genome_path, Chromosome_name, Reg_start, Reg_end, Edit_chrom_path)
