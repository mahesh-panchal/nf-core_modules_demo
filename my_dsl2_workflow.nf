#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// A workflow parameter called `reads`, which provides the path to a pair of Illumina sequence files.
params.reads = ''

// Include the FASTQC process definition from the module file.
include { FASTQC } from './modules/nf-core/software/fastqc' addParams(options:[:])

workflow {

    main:
    Channel.fromFilePairs(params.reads)
        .set { input_ch }

}
