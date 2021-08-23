#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Include the FASTQC process definition from the module file.
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:modules['fastqc'])

// Helper function to provide channel input in the correct format for nf-core modules.
def get_sample_info(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}

workflow {

    main:
    Channel.fromPath(params.input)
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { input_ch }

    FASTQC(input_ch)
    FASTQC.out.version.view()

}
