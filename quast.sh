#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 4
#SBATCH --mem  128GB


module load quast/5.0.2


source assembly/config.txt






for sample in "${input_list[@]}"; do
echo $sample &&
quast.py  assembly_output/pilon_assembly_nondedup/$sample\_bowtie_pilon_assembly.fasta \
     assembly_output_nexteratrim_clumped_pilon-assembly/$sample\_bwa_3trimmed_q20_clumped_pilon_assembly.fasta \
    -r ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta \
    -g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Orf50.gff \
    -1 assembly_fastq/$sample\1_3trimmed_q20_clumped.fastq.gz -2 assembly_fastq/$sample\2_3trimmed_q20_clumped.fastq.gz \
    -o assembly_output_nexteratrim_clumped/quast_$sample\output_pilon_assembly_bwa \
; done

# for sample in "${input_list[@]}"; do
# echo $sample &&
# quast.py assembly_output/pilon/$sample\_pilon_mh_1.fasta \
#     assembly_output/pilon/$sample\_pilon_sp_1.fasta  \
#     -r ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta \
#     -g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Orf50.gff \
#     -1 assembly_fastq/$sample\1_3trimmed_q20.fastq.gz -2 assembly_fastq/$sample\2_3trimmed_q20.fastq.gz \
#     -o assembly_output/quast_$sample\_pilon_1_output_2 \
# ; done




# for sample in "${input_list[@]}"; do
#     echo $sample &&
#     quast.py \
#         assembly_output/$sample\_Megahit_readassembly/final.contigs.fa \
#         assembly_output/pilon/$sample\_pilon_mh_1.fasta \
#         assembly_output/pilon/$sample\_pilon_mh_2.fasta \
#         assembly_output/pilon/$sample\_pilon_mh_3.fasta \
#         assembly_output/pilon/$sample\_pilon_mh_4.fasta \
#         assembly_output/pilon/$sample\_pilon_mh_5.fasta \
#         -r ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta \
#         -g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Orf50.gff \
#         -1 assembly_fastq/$sample\1_3trimmed_q20.fastq.gz -2 assembly_fastq/$sample\2_3trimmed_q20.fastq.gz \
#         -o assembly_output/quast_$sample\_pilon_iter_mh_output \
# ; done

# for sample in "${input_list[@]}"; do
#     echo $sample &&
#     quast.py \
#         assembly_output/$sample\_Spades_readassembly/scaffolds.fasta \
#         assembly_output/pilon/$sample\_pilon_sp_1.fasta \
#         assembly_output/pilon/$sample\_pilon_sp_2.fasta \
#         assembly_output/pilon/$sample\_pilon_sp_3.fasta \
#         assembly_output/pilon/$sample\_pilon_sp_4.fasta \
#         assembly_output/pilon/$sample\_pilon_sp_5.fasta \
#         -r ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta \
#         -g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Orf50.gff \
#         -1 assembly_fastq/$sample\1_3trimmed_q20.fastq.gz -2 assembly_fastq/$sample\2_3trimmed_q20.fastq.gz \
#         -o assembly_output/quast_$sample\_pilon_iter_sp_output \
# ; done