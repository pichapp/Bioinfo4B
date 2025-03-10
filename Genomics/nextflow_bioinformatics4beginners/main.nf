nextflow.enable.dsl=2

// Trimming and Filtering
// The large variation in quality scores should be handled before analysis.
// Step 1: Remove low-quality bases at the beginning and end of reads.
// Step 2: Filter reads with a high percentage of low-quality base calls.
// Parameters:
// - First 20 bases removed (-f 20)
// - Last position kept at 240 (-l 240)
// - Filter reads with at least 95% of bases above quality 30 (-q 30)
process TRIM_QUALITY_FILTER {
    publishDir params.outdir, mode: 'copy'  // Output of this process will be save to outdir
    input:
        path reads
    output:
        path "reads_filtered.fastq", emit: filtered_reads
        path "reads_trimmed.fastq", emit: trimmed_reads
    
    script:
    """
    fastx_trimmer -f 20 -l 240 -i ${reads} -o reads_trimmed.fastq
    fastq_quality_filter -q 30 -p 95 -i reads_trimmed.fastq -o reads_filtered.fastq
    """
}

// Quality Control (QC) Analysis
// Runs FastQC on raw, trimmed, and filtered reads to assess quality before and after processing.
process FASTQC_ANALYSIS {
    publishDir params.outdir, mode: 'copy'  // Output will be saved to the specified directory
    input:
        path reads
        path trimmed_reads
        path filtered_reads

    output:
        path "${reads.baseName}_fastqc", emit: raw_qc
        path "${trimmed_reads.baseName}_fastqc", emit: trimmed_qc
        path "${filtered_reads.baseName}_fastqc", emit: filtered_qc

        path "${reads.baseName}_fastqc.html", emit: raw_qc_html
        path "${trimmed_reads.baseName}_fastqc.html", emit: trimmed_qc_html
        path "${filtered_reads.baseName}_fastqc.html", emit: filtered_qc_html

    script:
    """
    fastqc ${reads}
    fastqc ${trimmed_reads}
    fastqc ${filtered_reads}
    
    unzip \$(basename ${reads} .fastq)_fastqc.zip
    unzip \$(basename ${trimmed_reads} .fastq)_fastqc.zip
    unzip \$(basename ${filtered_reads} .fastq)_fastqc.zip
    """
}


// Alignment to the Reference Genome
// Step 1: Index the reference genome using Bowtie2.
// Step 2: Align reads to the reference.
// Step 3: Convert the SAM file to BAM and sort it.
// Step 4: Index the sorted BAM file.
process ALIGNMENT {
    publishDir params.outdir, mode: 'copy'  // Output of this process will be save to outdir
    input:
        path reference
        path reads
    output:
        path "alignment_sorted.bam", emit: sorted_bam
        path "alignment.txt", emit: txt
    
    script:
    """
    bowtie2-build ${reference} ref
    bowtie2 -x ref -q ${reads} -S alignment.sam
    samtools faidx ${reference}
    samtools view -bt ${reference}.fai alignment.sam > alignment.bam
    samtools sort alignment.bam -o alignment_sorted.bam
    samtools index alignment_sorted.bam
    samtools tview -d T alignment_sorted.bam > alignment.txt
    """
}

process VISUALIZATION {
    publishDir params.outdir, mode: 'copy'  // Output of this process will be save to outdir
    input:
        path bam
        path reference
    output:
        path "depth.csv", emit: depth_csv
    
    script:
    """
    samtools depth ${bam} > depth.csv
    """
}

// Variant Calling
// Step 1: Generate pileup data with `bcftools mpileup`.
// Step 2: Call variants with `bcftools call`.
// Step 3: Generate variant statistics with `bcftools stats`.
process VARIANT_CALLING {
    publishDir params.outdir, mode: 'copy'  // Output of this process will be save to outdir
    input:
        path bam
        path reference
    output:
        path "calls.vcf"
        path "calls.vchk"
    
    script:
    """
    bcftools mpileup -f ${reference} ${bam} | bcftools call -mv -Ov -o calls.vcf
    bcftools stats calls.vcf > calls.vchk
    """
}

workflow {
    TRIM_QUALITY_FILTER(params.input)
    FASTQC_ANALYSIS(params.input, TRIM_QUALITY_FILTER.out.trimmed_reads, TRIM_QUALITY_FILTER.out.filtered_reads)
    ALIGNMENT(params.reference, TRIM_QUALITY_FILTER.out.filtered_reads)
    VISUALIZATION(ALIGNMENT.out.sorted_bam, params.reference)
    VARIANT_CALLING(ALIGNMENT.out.sorted_bam, params.reference)
}
