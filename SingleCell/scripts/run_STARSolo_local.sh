#!/usr/bin/env bash
set -eou pipefail

if [ -z ${OUTPUTDIR+x} ]; then
	>&2 echo "OUTPUTDIR was not set, defaulting to current directory '$(pwd)'";
	OUTPUTDIR=$(pwd);
fi
>&2 echo "Output will be generated in '${OUTPUTDIR}'"

if [ -z ${STARDATADIR+x} ]; then
	>&2 echo "STARDATADIR was not set, defaulting to current directory '$(pwd)'";
	STARDATADIR=$(pwd);
fi
>&2 echo "Looking for STAR and data files in '${STARDATADIR}'"

if [ -z ${NTHREADS+x} ]; then
	>&2 echo "NTHREADS was not set, defaulting to 4";
	NTHREADS=4;
fi

curl -L https://github.com/acg-team/Bioinfo4B/blob/main/SingleCell/data/STAR_data_archive.tar.gz?raw=true > ${STARDATADIR}/STAR_data_archive.tar.gz
tar -xzvf ${STARDATADIR}/STAR_data_archive.tar.gz -C ${STARDATADIR}

# Make a directory to store STAR output in
if [ ! -d ${OUTPUTDIR}/starsolo_out ]; then
	mkdir ${OUTPUTDIR}/starsolo_out
fi

# Run STAR in STARSolo mode
${STARDATADIR}/bin/STAR \
	--genomeDir ${STARDATADIR}/genome_idx/ \
	--readFilesCommand zcat \
	--readFilesIn ${STARDATADIR}/reads/subset100k_pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
				  ${STARDATADIR}/reads/subset100k_pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
	--runThreadN ${NTHREADS} \
	--outFileNamePrefix ${OUTPUTDIR}/starsolo_out/1kpmbc_ \
	--soloType Droplet \
	--soloCBwhitelist ${STARDATADIR}/barcodes/3M-february-2018.txt \
	--soloCBlen 16 \
	--soloUMIstart 17 \
	--soloUMIlen 12 
