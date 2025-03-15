#!/usr/bin/env bash
set -eou pipefail

if [ -z ${STARDATADIR+x} ]; then
	>&2 echo "STARDATADIR was not set, defaulting to current directory '$(pwd)'";
	STARDATADIR=$(pwd);
fi
>&2 echo "Looking for STAR in '${STARDATADIR}'"

if [ -z ${NTHREADS+x} ]; then
	>&2 echo "NTHREADS was not set, defaulting to 4";
	NTHREADS=4;
fi

# Get UCSC gene annotations for hg38
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz | gunzip > ${STARDATADIR}/hg38.knownGene.gtf

# Download UCSC assembly of hg38 chromosome 22
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | gunzip > ${STARDATADIR}/chr22.fa

# Generate index
mkdir ${STARDATADIR}/genome_idx
${STARDATADIR}/bin/STAR \
    --runThreadN ${NTHREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STARDATADIR}/genome_idx \
    --genomeFastaFiles ./chr22.fa \
    --sjdbGTFfile ./hg38.knownGene.gtf \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11