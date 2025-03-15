#!/usr/bin/env bash
#
# The #SBATCH lines below are //not// commented out! These lines are read by the Slurm preprocessor. They need to start with a '#' character to work.
#
#SBATCH --job-name=STARSolo
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00-00:30:00
#SBATCH --mem=15GB
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8

set -eou pipefail

if [ -z ${OUTPUTDIR+x} ]; then
	>&2 echo "OUTPUTDIR was not set, defaulting to home directory '${HOME}'";
	OUTPUTDIR=${HOME};
else
	if [ ! -d ${OUTPUTDIR} ]; then
		>&2 echo "The environment variable OUTPUTDIR must point to a directory where STRARSolo results can be written";
		exit 1;
	fi
fi

if [ -z ${STARDATADIR+x} ]; then
	STARDATADIR="/cfs/earth/scratch/shared/bioinfo4beginners/single_cell/";
fi
>&2 echo "Looking for STAR and data files in '${STARDATADIR}'"

# Change to our local scratch directory on the compute node, and make a directory to store STAR output in
cd ${LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH}
mkdir ./starsolo_out

# Run STAR in STARSolo mode
${STARDATADIR}/bin/STAR \
	--genomeDir ${STARDATADIR}/genome_idx/ \
	--readFilesCommand zcat \
	--readFilesIn ${STARDATADIR}/reads/subset100k_pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
		      ${STARDATADIR}/reads/subset100k_pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
	--runThreadN ${SLURM_CPUS_PER_TASK} \
	--outFileNamePrefix ./starsolo_out/1kpmbc_\
	--soloType Droplet \
	--soloCBwhitelist ${STARDATADIR}/barcodes/3M-february-2018.txt \
	--soloCBlen 16 \
	--soloUMIstart 17 \
        --soloUMIlen 12 


cp -r starsolo_out/1kpmbc_Solo.out ${OUTPUTDIR}


