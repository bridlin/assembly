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
output_dir=assembly_output_nexteratrim_clumped_pilon-assembly



assemblydic=assembly_output_nexteratrim_clumped_pilon-assembly
assemblyprefix=bwa_3trimmed_q20_clumped_pilon_assembly


echo $assemblydic
echo $output_dir
echo $fastq_directory

## 1. pilon round

mkdir -p $assemblydic\/bwa
mkdir -p $assemblydic\/pilon

for sample in "${input_list[@]}"; do
echo $sample &&
assembly=$assemblydic\/$sample\_$assemblyprefix\.fasta &&
echo $assembly &&
bwa-mem2 index -p  $assembly $assembly \
bwa-mem2 \
    mem $assembly $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq $fastq_directory/$sample\2_3trimmed_q20_clumped.fastq   \
    > $assemblydic\/bwa/$sample\aln-pe_$assemblyprefix\.sam &&
samtools view -S -b $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\.sam > $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\.sam.bam &&
samtools sort $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\.sam.bam -o $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\_sorted.bam &&
samtools index $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\_sorted.bam &&
rm -f  $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\.sam &&
rm -f  $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\.sam.bam &&
pilon \
    -Xmx8g \
    --genome $assembly \
    --bam $assemblydic\/bwa/$sample\aln-pe\_$assemblyprefix\_sorted.bam \
    --output $sample\_pilon_$assemblyprefix\_1 \
    --outdir $assemblydic\/pilon \
    --threads 4  \
    --changes \
    --tracks \
; done

# for sample in "${input_list[@]}"; do
# echo $sample &&
# assembly=$assemblydic\/$sample\_$assemblyprefix\.fasta &&
# echo $assembly &&
# bowtie2-build \
#     -f $assembly $assembly
# bowtie2 \
#     -x $assembly \
#     -1 $fastq_directory/$sample\1_3trimmed_q20_clumped.fastq \
#     -2 $fastq_directory/$sample\2_3trimmed_q20_clumped.fastq   \
#     -S $output_dir/bowtie2/$sample\aln-pe_$assemblyprefix\.sam \
#     2> $output_dir/bowtie2/$sample\_bowtie_$assemblyprefix\.log &&
# samtools view -S -b $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\.sam > $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\.sam.bam &&
# samtools sort $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\.sam.bam -o $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\_sorted.bam &&
# samtools index $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\_sorted.bam &&
# rm -f  $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\.sam &&
# rm -f  $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\.sam.bam &&
# pilon \
#     -Xmx8g \
#     --genome $assembly \
#     --bam $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\_sorted.bam \
#     --output $sample\_pilon_$assemblyprefix\_1 \
#     --outdir $output_dir/pilon \
#     --threads 4  \
#     --changes \
#     --tracks \
# ; done

#

### iterating pilon rounds

# for iter in "${iter_list[@]}"; do
#     y=$(($iter+1)) 
    
#     for sample in "${input_list[@]}"; do
#        assembly=$output_dir_pilon\/$sample\_pilon_$assemblyprefix\_$iter.fasta &&
#         echo $assembly &&
#         bowtie2-build \
#             -f $assembly $assembly
#         bowtie2 \
#             -x $assembly \
#             -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
#             -2 $fastq_directory/$sample\2_3trimmed_q20.fastq.gz   \
#             -S $output_dir/bowtie2/$sample\aln-pe_$assemblyprefix\_$iter\.sam \
#             2> $output_dir/bowtie2/$sample\_bowtie_$assemblyprefix\_$iter\.log &&
#         samtools view -S -b $output_dir/bowtie2/$sample\aln-pe_$assemblyprefix\_$iter\.sam > $output_dir/bowtie2/$sample\aln-pe_$assemblyprefix\_$iter\.sam.bam &&
#         samtools sort $output_dir/bowtie2/$sample\aln-pe_$assemblyprefix\_$iter\.sam.bam -o $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\_$iter\_sorted.bam
#         samtools index $output_dir/bowtie2/$sample\aln-pe\_$assemblyprefix\_$iter\_sorted.bam &&
#         rm -f  $output_dir_bowtie/$sample\aln-pe_$assemblyprefix\_$iter\.sam  &&
#         rm -f  $output_dir_bowtie/$sample\aln-pe_$assemblyprefix\_$iter\.sam.bam &&
#         pilon \
#             -Xmx8g \
#             --genome $assembly \
#             --bam $output_dir_bowtie/$sample\aln-pe\_$assemblyprefix\_$iter\_sorted.bam \
#             --output $sample\_pilon_$assemblyprefix\_$y \
#             --outdir $output_dir_pilon \
#             --threads 4  \
#             --changes \
#     ;done ;done    
        