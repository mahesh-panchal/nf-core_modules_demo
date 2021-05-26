#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// A workflow parameter called `reads`, which provides the path to a pair of Illumina sequence files.
params.reads = ''

workflow {

    main:
    Channel.fromFilePairs(params.reads)
        .set { input_ch }

}
