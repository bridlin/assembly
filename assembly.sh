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


module load cutadapt/4.0
module load trimmomatic/0.39
module load fastqc/0.11.9
module load megahit/1.2.9  
module load velvet/1.2.10
module load spades/3.15.5
module load multiqc/1.13

source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=$run\_output

mkdir $output_dir


for sample in "${input_list[@]}"; do
echo $sample &&
fastqc $fastq_directory/$sample\1.fastq.gz \
    --outdir $output_dir &&
fastqc $fastq_directory/$sample\2.fastq.gz \
    --outdir $output_dir &&
cutadapt  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
    -o $fastq_directory/$sample\1_3trimmed.fastq.gz \
    -p $fastq_directory/$sample\2_3trimmed.fastq.gz  \
    $fastq_directory/$sample\1.fastq.gz  $fastq_directory/$sample\2.fastq.gz \
    --minimum-length 40 \
    > $output_dir/$sample\_cutadapt_report.txt &&
trimmomatic PE \
    -threads 4 \
    -trimlog $output_dir/$sample\trim \
    $fastq_directory/$sample\1_3trimmed.fastq.gz $fastq_directory/$sample\2_3trimmed.fastq.gz \
    $fastq_directory/$sample\1_3trimmed_q20.fastq.gz   $fastq_directory/$sample\1_3trimmed_q20_un.fastq.gz $fastq_directory/$sample\2_3trimmed_q20.fastq.gz  $fastq_directory/$sample\2_3trimmed_q20_un.fastq.gz \
    SLIDINGWINDOW:4:20 \
    MINLEN:40 &&
fastqc $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
    --outdir $output_dir &&
fastqc $fastq_directory/$sample\2_3trimmed_q20.fastq.gz \
    --outdir $output_dir &&
megahit -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz -2 $fastq_directory/$sample\2_3trimmed_q20.fastq.gz -o $output_dir/$sample\_Megahit_readassembly  \
spades.py -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz -2 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz -o $output_dir/$sample\_Megahit_readassembly \
multiqc   \
    $output_dir \
    --outdir $output_dir 
;done