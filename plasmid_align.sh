#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 64GB


module load bowtie2/2.5.1
module load  samtools/1.13

source assembly/config.txt

fastq_directory=$run\_fastq
genome=plasmids/23pTB007_Cas9_T7_Tub.fa
output=plasmid_mapping

mkdir -p $output

for sample in "${input_list[@]}"; do
bowtie2-build \
    -f $genome $genome &&
bowtie2 -x $genome -f -p 8  -1 $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz  -2 $fastq_directory/$sample\2_3trimmed_q20_clumped.fastq.gz  -S $output\/$x\_plasmid_aligned.sam &&
samtools view -S -b  $output\/$x\_plasmid_aligned.sam  >  $output\/$x\_plasmid_aligned.bam  &&
samtools sort  $output\/$x\_plasmid_aligned.bam  -o  $output\/$x\_plasmid_aligned_sorted.bam  &&
samtools $output\/$x\_plasmid_aligned_sorted.bam   &&
rm -f $output\/$x\_plasmid_aligned.sam  &&
rm -f $output\/$x\_plasmid_aligned.bam  ; done

