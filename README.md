# nf-core Modules Demonstration

A demonstration of how to use nf-core modules in your workflow.
This repository include sample code as well as explanation on how
to use nf-core modules.

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
can be seen by following their respective links. nf-core provides a software package `nf-core`
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
is simply a collection of three files, `functions.nf`, `main.nf`, and `meta.yml`.
`functions.nf` contains utility functions for `main.nf`. The `main.nf` script
contains the process definition of the module. Finally, the `meta.yml` file
provides software metadata which can be used for input and output validation
among other things.

nf-core modules are written to a standard that makes it simpler to include them
across nf-core DSL2 workflows. However this makes it slightly more difficult to
include in a generic Nextflow DSL2 workflow.

To see how to use the nf-core modules in your workflow, follow the links below.

- [Minimal code example](./docs/MinimalCodeExample.md) - When code is easier to follow. More in depth explanation below.
- [Using an nf-core module](./docs/UsingNFCoreModules.md)
- [Adapting an nf-core module](./docs/AdaptingNFCoreModules.md)
- [Making a new local module](./docs/MakingALocalModule.md)
- [Advanced module parameterisation](./docs/AdvancedParameterisation.md)
