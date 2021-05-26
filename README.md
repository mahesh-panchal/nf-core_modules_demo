# nf-core Modules Demonstration

A demonstration of how to use nf-core modules in your workflow.

[nf-core](https://nf-co.re/) is a set of community curated workflows
written in Nextflow. The DSL2 syntax of Nextflow allows the use of
Modules, scripts that can be included and shared across workflows.

For example,

```nextflow
include { foo } from './some/module'
```
where `foo` is a `process` in the `module.nf` script.

nf-core introduced a [modules repository](https://github.com/nf-core/modules/)
to aid the construction of workflows using the DSL2 syntax. These
nf-core modules are intended to provide software tool specific processes
which can be reused across nf-core workflows. For example,
[FastQC](https://github.com/nf-core/modules/tree/master/software/fastqc),
[BWA index](https://github.com/nf-core/modules/tree/master/software/bwa/index), and
[BWA mem](https://github.com/nf-core/modules/tree/master/software/bwa/mem) implementations
can be seen by following their respective links. nf-core provides a software `nf-core`
(via `conda` or `pip`) that provides tools to aid using nf-core modules, among other things.

Modules available in the nf-core modules repository can be listed using
```bash
nf-core modules list
```

A module can be included using
```bash
nf-core modules install --tool <module> <directory>
```
For example,
```bash
nf-core modules install --tool fastqc .
```
which installs the `fastqc` module in the current directory `.`. The module
is simply a collection of three files, `functions.nf`, `main.nf`, and `meta.yml`,
explained later.

nf-core modules are written to a standard that makes it easy to include them
across nf-core DSL2 workflows. However this makes it slightly more difficult to
include in a generic Nextflow DSL2 workflow.

## How to use an nf-core module

Let's assume you've started a DSL2 workflow like so, called `my_dsl2_workflow.nf`:

```nextflow
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
