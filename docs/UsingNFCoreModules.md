## How to use an nf-core module

Let's assume you've started a DSL2 workflow like so, called `my_dsl2_workflow.nf`:

```nextflow
#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// A workflow parameter called `input`, which provides the path to CSV samplesheet.
params.input = ''

workflow {

    main:
    Channel.fromPath(params.input)
        .set { input_ch }

}
```

You've discovered using `nf-core modules list` that software you want, FastQC, is available as an
nf-core module. So let's install it in the same directory as your DSL2 workflow.

```bash
nf-core modules install --tool fastqc .
```

Using this the first time around will likely produce an error, along the lines of
`CRITICAL Could not find a 'main.nf' or 'nextflow.config' file in '.'`. To remedy
this, create a `nextflow.config` file (`touch nextflow.config`), and then try to
install the module again.

Your working directory will now look something like this:
```
| - modules/nf-core/software/fastqc
|    | - functions.nf
|    | - main.nf
|    \ - meta.yml
|
| - my_dsl2_workflow.nf
\ - nextflow.config
```

Let's now include the `FASTQC` process from the fastqc module to our DSL2 workflow.

```nextflow
#! /usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl = 2

// A workflow parameter called `input`, which provides the path to CSV samplesheet.
params.input = ''

// Include the FASTQC process definition from the module file.
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:[:])

workflow {

    main:
    Channel.fromPath(params.input)
        .set { input_ch }

}
```

Let's unpack the line we just included.

```nextflow
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:[:])
```

An nf-core module convention is to name processes in uppercase using `<SOFTWARE>_<TOOL>`.
This means our FastQC process will be named `FASTQC`, which we can see in the
`./modules/nf-core/software/fastqc/main.nf` file. The path following the `from` keyword
tells Nextflow to look for the process definition in the path `./modules/nf-core/software/fastqc/main`.
Nextflow appends `.nf` to the path and checks the path to see if the file `main.nf` exists.
The last part is the `addParams`
call, which passes parameters to the workflow (`main.nf`) that contains the process
definition. In this case, `addParams(options:[:])` initialises the `main.nf` workflow parameter
`params.options` to an empty Map `[:]`.

### Providing module parameters.

Many nf-core modules may need additional parameters.
As we've just read, parameters are passed to the module via the `addParams(Map params)` call.
The convention for nf-core Modules and DSL2 workflows is to provide module parameters in a
file called `conf/modules.config`, which is included by the `nextflow.config`. Let's add
those as a first step.

```bash
mkdir conf
touch conf/modules.config
```
and at the top of the `nextflow.config` add:
```nextflow
include 'conf/modules.config'
```

The nf-core DSL2 template `conf/modules.config` is formatted in a particular way.
There is a large comment at the top of the file describing which parameters can
be passed and what they do. This is followed by a `params` block, which includes
a `modules` block, and can be accessed in the DSL2 workflow using `params.modules`.

```
/*
========================================================================================
    Config file for defining DSL2 per module options
========================================================================================
    Available keys to override module options:
        args            = Additional arguments appended to command in module.
        args2           = Second set of arguments appended to command in module (multi-tool modules).
        args3           = Third set of arguments appended to command in module (multi-tool modules).
        publish_dir     = Directory to publish results.
        publish_by_meta = Groovy list of keys available in meta map to append as directories to "publish_dir" path
                            If publish_by_meta = true                 - Value of ${meta['id']} is appended as a directory to "publish_dir" path
                            If publish_by_meta = ['id', 'custompath'] - If "id" is in meta map and "custompath" isn't then "${meta['id']}/custompath/"
                                                                        is appended as a directory to "publish_dir" path
                            If publish_by_meta = false / null         - No directories are appended to "publish_dir" path
        publish_files   = Groovy map where key = "file_ext" and value = "directory" to publish results for that file extension
                            The value of "directory" is appended to the standard "publish_dir" path as defined above.
                            If publish_files = null (unspecified)     - All files are published.
                            If publish_files = false                  - No files are published.
        suffix          = File name suffix for output files.
----------------------------------------------------------------------------------------
*/

params {
    modules {
    }
}
```

To add parameters for a tool, a named block is included in the `modules` block with the parameters to include/override.
```
params {
    modules {
        'fastqc' {
            args         = '--quiet'
            publishDir   = '01_FastQC'
        }
    }
}
```
The convention for naming the module options block is generally `<subworkflow>_<software>_<tool>`,
however any string can be used.

The nf-core DSL2 template includes an extra line to help accessing these parameter blocks.

```nextflow
// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
```

By adding the line above to your DSL2 workflow, you can then pass the `fastqc` module options to the FastQC module, like so:
```
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:modules['fastqc'])
```
This passes the Map `modules['fastqc']` as the `params.options` of the FastQC `main.nf` script.
You can see where the parameters are used in the
`main.nf` scripts by looking for `$options.<parameter>`, for example `$options.args`.

Amend your workflow to include these changes:
```nextflow
// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Include the FASTQC process definition from the module file.
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams(options:modules['fastqc'])
```

### Providing metadata.

If you looked at the `main.nf` script of the FastQC module, you
may have noticed that the input declaration looks like:
```
    input:
    tuple val(meta), path(reads)
```
This means the `FASTQC` process is looking for input in the form
`[meta, reads]` where `meta` is sample metadata, and `reads` are
the read files to be processed. `meta` is a Map which provides
sample specific information, such as `id`, `single_end`, and perhaps other
fields e.g. `read_group`. Each module needs to be individually checked for
which sample metadata fields need to be provided (search `main.nf` for
variables named `meta.<field>`, e.g. `meta.id`). Not all nf-core
modules require sample metadata, for example
[BWA INDEX](https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf).

In the nf-core DSL2 template, the sample metadata is expected to be set by
the workflow `INPUT_CHECK`, which calls the process `SAMPLESHEET_CHECK`,
which runs the python script `check_samplesheet.py` from the workflow `bin` folder.
The nf-core DSL2 template `check_samplesheet.py` modifies the input sample sheet
to contain the columns `sample`,`single_end`,`fastq_1`,`fastq_2`, after being
provided `sample`,`fastq_1`, and `fastq_2` in the samplesheet. The `INPUT_CHECK`
workflow then uses the `splitCsv` and `map` channel operators to construct the sample metadata
using the `get_sample_info` function, which returns the complex data structure
`[{id: <string>, single_end: <boolean>}, [fastq_1, fastq_2]]`. This complex data
structure is an Array with two elements. The first element is a Map with the keys
`id`, and `single_end`. The second element is an Array containing the paths
to the sample input files.

If you're not using the nf-core DSL2 template, you can use a similar method in your workflow to
provide sample metadata. Provide a samplesheet with the sample metadata as extra
columns, and then using the `splitCsv`, and `map` channel operators to convert the
input to the correct format.

For example, if the samplesheet contains the columns `sample`,`single_end`,`fastq_1`,
and `fastq_2`, you could define your own function based on the `get_sample_info`
function from `INPUT_CHECK` to provide channel data in the expected format.
If you're comfortable with Nextflow DSL2 syntax, you could squirrel the function away in a
utility file for tidiness.

```nextflow
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
```

and then in the workflow do:
```nextflow
workflow {

    main:
    Channel.fromPath(params.input)
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { input_ch }

}
```

Your channel data should now be in a format nf-core modules expect, allowing you to use
the modules in your workflow.

```nextflow
workflow {

    main:
    Channel.fromPath(params.input)
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { input_ch }

    FASTQC(input_ch)

}
```

### Providing parameter defaults for modules

Before the workflow will function, the nf-core modules still
require further parameterisation. The `main.nf` script for the
FastQC module still expects values for the variables
`params.enable_conda`, `params.outdir`, `params.publish_dir_mode`,
and `params.singularity_pull_docker_container`. You
can do this by adding the following block to your `nextflow.config`.

```nextflow
params {
    enable_conda = false
    outdir = './results'
    publish_dir_mode = 'copy'
    singularity_pull_docker_container = true
}
```

This `params` block provides default values for the variables above.
`params` is treated as a global variable, which is why it is available
to the module scripts without passing them through the `addParams` call.

### Testing the workflow

The workflow should now be able to use the modules and complete successfully.
You can test this by creating a `test_samplesheet.csv`.

```
sample,single_end,fastq_1,fastq_2
test,0,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
```

```bash
nextflow run my_dsl2_workflow.nf --input test_samplesheet.csv
```
