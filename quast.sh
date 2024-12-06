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
./quast.py assembly_output/$sample\_Megahit_readassembly/final.contigs.fa \
            assembly_output/$sample\_Spades_readassembly/scaffolds.fasta \
            assembly_output/$sample\_Spades_readassembly/contigs.fasta \
        -r ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LmexicanaMHOMGT2001U1103/TriTrypDB-68_LmexicanaMHOMGT2001U1103_Genome.fasta \
        -g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LmexicanaMHOMGT2001U1103/TriTrypDB-68_LmexicanaMHOMGT2001U1103_Orf50.gff \
        -1 assembly_fastq/$sample\1_3trimmed_q20.fastq.gz -2 assembly_fastq/$sample\2_3trimmed_q20.fastq.gz \
        -o assembly_output/quast_$sample\output
; done