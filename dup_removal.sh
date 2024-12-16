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
module load multiqc/1.13
module load bbmap/39.00
module load fastqc/0.11.9

source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=assembly_output_2

mkdir -p $output_dir

for sample in "${input_list[@]}"; do
    cutadapt  -a CTGTCTCTTATACACATCT -a AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -A AGATGTGTATAAGAGACAG \
        -o $fastq_directory/$sample\1_3trimmed_2.fastq.gz \
        -p $fastq_directory/$sample\2_3trimmed_2.fastq.gz  \
        -n 3 \
        $fastq_directory/$sample\1.fastq.gz  $fastq_directory/$sample\2.fastq.gz \
        --minimum-length 40 \
        > $output_dir/$sample\_cutadapt_report_2.txt &&
    fastqc $fastq_directory/$sample\1_3trimmed.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\1_3trimmed_2.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed_2.fastq.gz \
        --outdir $output_dir &&
    ### clumping reads
    clumpify.sh \
        in1=$fastq_directory/$sample\1_3trimmed.fastq.gz \
        in2=$fastq_directory/$sample\2_3trimmed.fastq.gz \
        out1=$fastq_directory/$sample\1_3trimmed_clumped.fastq.gz \
        out2=$fastq_directory/$sample\2_3trimmed_clumped.fastq.gz \
        dedupe=t \
        optical=f &&
    clumpify.sh \
        in1=$fastq_directory/$sample\1_3trimmed_2.fastq.gz \
        in2=$fastq_directory/$sample\2_3trimmed_2.fastq.gz \
        out1=$fastq_directory/$sample\1_3trimmed_2_clumped.fastq.gz \
        out2=$fastq_directory/$sample\2_3trimmed_2_clumped.fastq.gz \
        dedupe=t \
        optical=f &&
    fastqc $fastq_directory/$sample\1_3trimmed_clumped.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed_clumped.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\1_3trimmed_2_clumped.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed_2_clumped.fastq.gz \
        --outdir $output_dir \
;done

multiqc   \
    $output_dir \
    --outdir $output_dir  \