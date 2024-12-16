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

source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=assembly_output_2

mkdir -p $output_dir

for sample in "${input_list[@]}"; do
    cutadapt  -a CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA -A CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA \
        -o $fastq_directory/$sample\1_3trimmed.fastq.gz \
        -p $fastq_directory/$sample\2_3trimmed.fastq.gz  \
        $fastq_directory/$sample\1.fastq.gz  $fastq_directory/$sample\2.fastq.gz \
        --minimum-length 40 \
        > $output_dir/$sample\_cutadapt_report.txt &&
    fastqc $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed_q20.fastq.gz \
        --outdir $output_dir &&
    ### clumping reads
    clumpify.sh \
        in1=$fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
        in2=$fastq_directory/$sample\2_3trimmed_q20.fastq.gz \
        out1=$fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz \
        out2=$fastq_directory/$sample\2_3trimmed_q20_clumped.fastq.gz \
        dedupe=t \
        optical=f \
        outd=$fastq_directory/$sample\duplicates.fq &&
    fastqc $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz \
        --outdir $output_dir &&
    fastqc $fastq_directory/$sample\2_3trimmed_q20_clumped.fastq.gz \
        --outdir $output_dir 
;done

multiqc   \
    $output_dir \
    --outdir $output_dir  \