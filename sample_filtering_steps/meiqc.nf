#!/usr/bin/env nextflow


params.samplesheet = "$baseDir/samplesheet_example.csv"
// samplesheet uses germline SNV/indel vcfs

params.fasta_file = "$baseDir/hg38.fa"
params.fasta_file_idx = "$baseDir/hg38.fa.fai"

params.aeon_afs= "<path>/aeon/refs/g1k_allele_freqs.txt"
params.aeon_labs = "<path>/aeon/refs/pop2super.txt"

params.outdir = "$baseDir/results"

workflow {

    sample_channel = Channel
                        .fromPath( params.samplesheet )
                        .splitCsv( header: true, sep: ',' )
                        .map { row -> tuple( row.sample_id, file(row.bam_file), file(row.bai_file), file(row.vcf_file), file(row.tbi_file) ) }

    bam_ch = channel.fromPath(params.bam_files)
            .splitText(){ it.trim() }
            .map { file(it) }

    bai_ch = channel.fromPath(params.bai_files)
            .splitText(){ it.trim() }
            .map { file(it) }

    buffered_bam_ch = bam_ch.buffer( size: 10, remainder: true )

    PANDEPTH(sample_channel,params.outdir)
    INDEXCOV(sample_channel,params.fasta_file_idx,params.outdir)
    ANCESTRY(sample_channel,params.aeon_afs,params.aeon_labs,params.outdir)
}

process ANCESTRY {
    tag "aeon_$sample_id"

    cpus 1
    memory '1000M'
    time '60m'

    publishDir "${params.outdir}/aeon", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(vcf_file), path(tbi_file)
    path afs
    path labs
    path outdir

    output:
    path "${sample_id}_pop.txt"

    script:
    """
    aeon.py $vcf_file -o ${sample_id} --no_visualisation --population_labels $labs --allele_freqs $afs
    parse_aeon_multi_sample.py ${sample_id}_ae.csv > ${sample_id}_pop.txt
    """
}


process PANDEPTH {
    tag "pandepth_$sample_id"

    cpus 3
    memory '2G'
    time '30m'

    publishDir "${params.outdir}/pandepth", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(vcf_file), path(tbi_file)
    path outdir

    output:
    path "${sample_id}_mean_cov.txt"

    script:
    """
    pandepth -i ${bam_file} -o ${sample_id}
    zcat ${sample_id}.chr.stat.gz | tail -n 1 | awk '{print \$8}' - > ${sample_id}_mean_cov.txt
    """
}

process INDEXCOV {
    tag "indexcov_$sample_id"

    cpus 1
    memory '500MB'
    time '30m'

    publishDir "${params.outdir}/indexcov", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(vcf_file), path(tbi_file)
    path fasta_file_idx
    path outdir

    output:
    path "${sample_id}/${sample_id}.indexcov.tsv"
    path "${sample_id}/${sample_id}-indexcov.ped"

    script:
    """
    goleft indexcov --extranormalize -d $sample_id -fai ${fasta_file_idx} $bai_file
    parse_indexcovbed.py ${sample_id}/${sample_id}-indexcov.bed.gz ${sample_id}/${sample_id}.indexcov.tsv
    """
}

