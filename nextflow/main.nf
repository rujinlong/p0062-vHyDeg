#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """\
NF - vHyDeg PIPELINE
=========================
result: ${params.outdir}
report: ${params.report}
"""

include { hmmscan_imgvr } from 'modules/annotation'

workflow vhydeg {
    ch_gene = Channel.fromPath("${params.input_contigs}", checkIfExists: true).first()

    // MODULE: CheckV
    hmmscan_imgvr(ch_gene)

    // ANNOTATION (AMG)
    if ( params.use_dram ) {
        DRAMV (ch_vs2contigs, VIRSORTER2.out.vs2_affi_ch)
    }
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
