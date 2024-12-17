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

echo "test"
module load bowtie2/2.5.1
module load samtools/1.21
module load pilon/1.23
module load bwa-mem2/2.2.1

echo "test2"

source assembly/config.txt

echo "test3"

fastq_directory=$run\_fastq
output_dir=$run\_output_nexteratrim_clumped
genome='../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta'

echo "test4"

output_bwa=$output_dir\_bwa 
output_bowtie=$output_dir\_bowtie2
output_pilon=$output_dir\_pilon-assembly

mkdir -p $output_bwa
mkdir -p $output_bowtie
mkdir -p $output_pilon

echo "test5"


echo $genome

# bowtie2-build -f $genome $genome 

for sample in "${input_list[@]}"; do
echo "test6"
echo $sample &&
bowtie2 \
    -x $genome \
    -1 $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz \
    -2 $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz   \
    -S $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped.sam \
    2> $output_bowtie/$sample\_bowtie_3trimmed_q20_clumped.log &&
samtools view -S -b $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\.sam > $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam &&
samtools sort $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam -o $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam &&
samtools index $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam &&
rm -f  $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\.sam &&
rm -f  $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam &&
pilon \
    -Xmx32g \
    --genome $genome \
    --bam $output_bowtie/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam \
    --output $sample\_bowtie_3trimmed_q20_clumped_pilon_assembly \
    --outdir $output_pilon \
    --threads 4  \
    --changes \
    --tracks \
    --vcf \
    --fix snps,indels \
; done

# bwa-mem2 index  -p TriTrypDB-68_LinfantumJPCM5_Genome $genome

for sample in "${input_list[@]}"; do
echo $sample &&
bwa-mem2 \
    mem $genome  $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq.gz  $fastq_directory/$sample\2_3trimmed_q20_clumped.fastq.gz \
    > $output_bwa/$sample\aln-pe_3trimmed_q20_clumped.sam &&     
samtools view -S -b $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\.sam > $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam &&
samtools sort $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam -o $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam &&
samtools index $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam &&
rm -f  $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\.sam &&
rm -f  $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\.sam.bam &&
pilon \
    -Xmx100g \
    --genome $genome \
    --bam $output_bwa/$sample\aln-pe_3trimmed_q20_clumped\_sorted.bam \
    --output $sample\_bwa_3trimmed_q20_clumped_pilon_assembly \
    --outdir $output_pilon \
    --threads 4  \
    --changes \
    --tracks \
    --vcf \
    --fix snps,indels \
; done