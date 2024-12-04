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


source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=kraken2-results_$run\_5prime-trimmed

mkdir $output_dir


for sample in "${input_list[@]}"; do
echo $sample &&
fastqc $fastq_directory/$sample\L001_R1_001.fastq.gz \
    --outdir $output_dir &&
fastqc $fastq_directory/$sample\L001_R2_001.fastq.gz \
    --outdir $output_dir &&
cutadapt  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
    -o $fastq_directory/$sample\L001_R1_001_3trimmed.fastq.gz \
    -p $fastq_directory/$sample\L001_R2_001_3trimmed.fastq.gz  \
    $fastq_directory/$sample\L001_R1_001.fastq.gz  $fastq_directory/$sample\L001_R2_001.fastq.gz \
    --minimum-length 40 \
    > $output_dir/$sample\_cutadapt_report.txt &&
trimmomatic PE \
    -threads 4 \
    -trimlog $output_dir/$sample\trim \
    $fastq_directory/$sample\L001_R1_001_3trimmed.fastq.gz $fastq_directory/$sample\L001_R2_001_3trimmed.fastq.gz \
    $fastq_directory/$sample\L001_R1_001_3trimmed_q20.fastq.gz   $fastq_directory/$sample\L001_R1_001_3trimmed_q20_un.fastq.gz $fastq_directory/$sample\L001_R2_001_3trimmed_q20.fastq.gz  $fastq_directory/$sample\L001_R2_001_3trimmed_q20_un.fastq.gz \
    SLIDINGWINDOW:4:20 \
    MINLEN:40 &&
fastqc $fastq_directory/$sample\L001_R1_001_3trimmed_q20.fastq.gz \
    --outdir $output_dir &&
fastqc $fastq_directory/$sample\L001_R2_001_3trimmed_q20.fastq.gz \
    --outdir $output_dir &&