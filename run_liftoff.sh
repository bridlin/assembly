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


module load liftoff/1.6.3


liftoff \
-g ../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5.gff \
-chroms chroms_lift.txt   \
-o  LinS9F1T7Cas9.gff  \
-f features.txt \
assembly_output_nexteratrim_clumped_pilon-assembly_fix-all/pilon/LinS9F1T7Cas9.fasta \
../chu_diag_microbiom_setup/genome/TriTrypDB-68_LinfantumJPCM5/TriTrypDB-68_LinfantumJPCM5_Genome.fasta