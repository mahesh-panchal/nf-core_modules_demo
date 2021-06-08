## Minimal code example

The example here shows the minimum code necessary to use an nf-core
module in your workflow.

The workflow takes a formatted samplesheet and runs FastQC.

`main.nf`:
```nextflow
#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// A workflow parameter called `input`, which provides the path to CSV samplesheet.
params.input = ''

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Installed with `nf-core modules install -t fastqc .`
// Include the FASTQC process definition from the module file.
// modules['fastqc'] is Groovy Map containing parameters from the `modules.config`.
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:modules['fastqc'])

// Helper function to provide channel input in the correct format for nf-core modules.
def get_sample_info(LinkedHashMap row) {
    // This function is applied to every row of the input CSV.
    def meta = [:]                                        // Empty Groovy Map
    meta.id           = row.sample                        // Assign value of "sample" column to "id".
    meta.single_end   = row.single_end.toBoolean()        // Assign value of "single_end" to "single_end".

    def array = []                                        // Empty Groovy List
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        // Make a Groovy List with two elements.
        // First element is the Groovy Map "meta".
        // Second element is a Groovy List of Nextflow file objects.
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    // returns data entry in the format expected by FastQC module
    // FASTQC input declaration expects `tuple val(meta), path(reads)`
    // array is in the format `[ meta, reads ]`.
    return array
}

workflow {

    main:
    Channel.fromPath(params.input)
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { input_ch }

    FASTQC(input_ch)

}
```

`nextflow.config`:
```nextflow
// Load modules.config for DSL2 module specific options
includeConfig 'conf/modules.config'

// Include minimum necessary parameters for nf-core modules
params {
    enable_conda = false
    outdir = './results'
    publish_dir_mode = 'copy'
    singularity_pull_docker_container = true
}
```

`conf/modules.config`:
```nextflow
params {
    modules {
        'fastqc' {       // passed as modules['fastqc'] in the module include statement
            args         = '--quiet'
            publishDir   = '01_FastQC'
        }
    }
}
```

`samplesheet.csv`:
```
sample,single_end,fastq_1,fastq_2
test,0,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
```
The column header names are the keywords used to access the row
values in the helper function `get_sample_info` in the `main.nf` script.
