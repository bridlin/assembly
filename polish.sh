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


module load bowtie2/2.5.1
module load samtools/1.21
module load pilon/1.23

source assembly/config.txt

fastq_directory=$run\_fastq
output_dir=$run\_output

output_dir_bowtie=$run\_output/bowtie2
output_dir_pilon=$run\_output/pilon
mkdir $output_dir_bowtie
mkdir $output_dir_pilon

megahit_assembly=$output_dir/F6__Megahit_readassembly/final.contigs.fa
spades_assembly=$output_dir/F6__Spades_readassembly/scaffolds.fasta

echo $megahit_assembly
echo $spades_assembly


for sample in "${input_list[@]}"; do
echo $sample &&
bowtie2-build \
    -f $megahit_assembly $megahit_assembly
bowtie2 \
    -x $megahit_assembly \
    -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
    -2 $fastq_directory/$sample\2_3trimmed_q20.fastq.gz   \
    -S $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam \
    2> $output_dir/bowtie/reports/$sample\_bowtie_mh.log &&
samtools view -S -b $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam > $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam.bam &&
samtools sort $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam.bam -o $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh_sorted.bam &&
samtools index $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh_sorted.bam &&
rm -f  $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam &&
rm -f  $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh.sam.bam &&
pilon --genome $megahit_assembly --bam $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_mh_sorted.bam --output $output_dir/pilon/$sample\_pilon_mh --outdir $output_dir/pilon --threads 4  --changes --tracks \
; done

for sample in "${input_list[@]}"; do
echo $sample &&
bowtie2-build \
    -f $spades_assembly $spades_assembly
bowtie2 \
    -x $spades_assembly \
    -1 $fastq_directory/$sample\1_3trimmed_q20.fastq.gz \
    -2 $fastq_directory/$sample\2_3trimmed_q20.fastq.gz   \
    -S $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam \
    2> $output_dir/bowtie/reports/$sample\_bowtie_sp.log &&
samtools view -S -b $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam > $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam.bam &&
samtools sort $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam.bam -o $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp_sorted.bam &&
samtools index $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp_sorted.bam &&
rm -f  $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam &&
rm -f  $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp.sam.bam &&
pilon --genome $spades_assembly --bam $output_dir/bowtie/$sample\aln-pe_$genome_prefix\_sp_sorted.bam --output $output_dir/pilon/$sample\_pilon_sp --outdir $output_dir/pilon --threads 4  --changes --tracks \
; done