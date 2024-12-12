#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user bridlin.barckmann@umontpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 4
#SBATCH --mem  128GB


module load bowtie2/2.5.1
module load samtools/1.21
module load pilon/1.23

source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=$run\_output
genome='../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta'

mkdir -p $output_dir/bowtie2
mkdir -p $output_dir/bwa
mkdir -p $output_dir/pilon_assembly

echo $genome

# bowtie2-build -f $genome $genome 

for sample in "${input_list[@]}"; do
echo $sample &&
# bowtie2 \
#     -x $genome \
#     -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
#     -2 $fastq_directory/$sample\2_3trimmed_q20.fastq.gz   \
#     -S $output_dir/bowtie2/$sample\aln-pe.sam \
#     2> $output_dir/bowtie2/$sample\_bowtie.log &&
# samtools view -S -b $output_dir/bowtie2/$sample\aln-pe\.sam > $output_dir/bowtie2/$sample\aln-pe\.sam.bam &&
# samtools sort $output_dir/bowtie2/$sample\aln-pe\.sam.bam -o $output_dir/bowtie2/$sample\aln-pe\_sorted.bam &&
# samtools index $output_dir/bowtie2/$sample\aln-pe\_sorted.bam &&
# rm -f  $output_dir/bowtie2/$sample\aln-pe\.sam &&
# rm -f  $output_dir/bowtie2/$sample\aln-pe\.sam.bam &&
pilon \
    -Xmx32g \
    --genome $genome \
    --bam $output_dir/bowtie2/$sample\aln-pe\_sorted.bam \
    --output $sample\_bowtie_pilon_assembly \
    --outdir $output_dir/pilon_assembly \
    --threads 4  \
    --changes \
    --tracks \
    --vcf \
    --fix snps,indels \
; done

# bwa index $genome 

# for sample in "${input_list[@]}"; do
# echo $sample &&
#  bwa-mem2 \
#     mem $genome \  
#     $fastq_directory/$sample\1_3trimmed_q20.fastq.gz\
#     $fastq_directory/$sample\2_3trimmed_q20.fastq.gz \
#     > $output_dir/bwa/$sample\aln-pe.sam &&     
# samtools view -S -b $output_dir/bwa/$sample\aln-pe\.sam > $output_dir/bwa/$sample\aln-pe\.sam.bam &&
# samtools sort $output_dir/bwa/$sample\aln-pe\.sam.bam -o $output_dir/bwa/$sample\aln-pe\_sorted.bam &&
# samtools index $output_dir/bwa/$sample\aln-pe\_sorted.bam &&
# rm -f  $output_dir/bwa/$sample\aln-pe\.sam &&
# rm -f  $output_dir/bwa/$sample\aln-pe\.sam.bam &&
# pilon \
#     -Xmx100g \
#     --genome $genome \
#     --bam $output_dir/bwa/$sample\aln-pe\_sorted.bam \
#     --output $sample\_bwa_pilon_assembly \
#     --outdir $output_dir/pilon_assembly \
#     --threads 4  \
#     --changes \
#     --tracks \
#     --vcf \
#     --fix snps,indels \
# ; done