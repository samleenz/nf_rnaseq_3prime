// a nextflow pipeline for processing 3' bulk rnaseq (lexogen)

// written by sam lee
// TPH, ONJCRI


// load modules
include { RUN_FASTQC } from "${params.module_dir}/fastqc"
include { RUN_MULTIQC } from "${params.module_dir}/multiqc"
include { RUN_TRIM_READS_SE } from "./modules/cutadapt-se" 
include { RUN_STAR_GENOME_GENERATE } from "${params.module_dir}/star-index" 
include { RUN_STAR_ALIGN_SE } from "./modules/star-align-se"
include { RUN_STAR_ALIGN_SE as RUN_STAR_ALIGN_SE_2PASS } from "./modules/star-align-se"


// workflow
workflow {

    reads_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map(row -> tuple(row.sample, [file(row.read1)] ))

    
    trimmed_reads_ch = RUN_TRIM_READS_SE(reads_ch)

    fastqc_ch = RUN_FASTQC(trimmed_reads_ch)



    // STAR ALIGNMENT -------

    // make star genome index
    star_idx_ch = RUN_STAR_GENOME_GENERATE(
        params.genome, // genome
        params.gtf, // gtf
        params.tx_ver, // ensembl version
        params.read_len // read_length
    )
  

    // pass one
    (star_align_pass1_ch, star_junctions_pass1_ch, star_logs_pass1_ch) = RUN_STAR_ALIGN(
        trimmed_reads_ch,
        star_idx_ch,
        params.read_len,
        file(params.juncs), // see https://nextflow-io.github.io/patterns/optional-input/ for description of logic
        1
    )
    

    star_junctions = star_junctions_pass1_ch.collect( sort : true )

    // star_junctions.view()

    // pass two
    (star_align_pass2_ch, star_junctions_pass2_ch, star_logs_pass2_ch) = RUN_STAR_ALIGN_2PASS(
        trimmed_reads_ch,
        star_idx_ch,
        params.read_len,
        star_junctions_pass1_ch.collect( sort: true ),
        2
    )


    // FEATURECOUNTS -------- 

    (fc_counts_ch, fc_logs_ch) = RUN_FEATURECOUNTS(
        star_align_pass2_ch,
        params.gtf
    )

    fc_combine_ch = RUN_FEATURECOUNTS_COMBINE(fc_counts_ch.collect())

    make_se_ch = RUN_SUMMARIZED_EXPERIMENT(
        fc_combine_ch//,
        // params.samplesheet,
        // params.gtf,
        // params.project
    )

    // MULTIQC --------------

    multiqc_ch = RUN_MULTIQC(
        fastqc_ch
            .concat(star_logs_pass2_ch)
            .concat(fc_logs_ch)
            .collect()
        )



}

workflow.onComplete {
    log.info( 
        workflow.success ? 
        "Pipeline Complete!\n" : 
        "Pipeline Failed.\n" )
}