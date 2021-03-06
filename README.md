# nf-core Modules Demonstration

Tested with:
- Nextflow: 21.04.0
- nf-core: 2.1

A demonstration of how to use nf-core modules in your workflow.
This repository includes sample code as well as an explanation on how
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
[FastQC](https://github.com/nf-core/modules/tree/master/modules/fastqc),
[BWA index](https://github.com/nf-core/modules/tree/master/modules/bwa/index), and
[BWA mem](https://github.com/nf-core/modules/tree/master/modules/bwa/mem) implementations
can be seen by following their respective links. nf-core provides a software package `nf-core`
(via `conda` or `pip`) that provides tools to aid using nf-core modules, among other things.

To start using an nf-core module, some conditions need to be satisfied.
1. The directory you're installing to must have a `main.nf` script.
1. The directory you're installing to must have a  `nextflow.config`.
1. The `nextflow.config` must contain a `manifest` block with the
`name` and `homePage` parameters set.
1. The directory you're installing to must have a `modules` directory.
1. `nextflow` must be available in the `PATH` of the current environment.

Modules available in the nf-core modules repository can be listed using
```bash
nf-core modules list
```

A module can be included using
```bash
nf-core modules install <module>
```
For example,
```bash
nf-core modules install fastqc
```
which installs the `fastqc` module in a subdirectory `./modules/nf-core/modules/fastqc/`. The module
is a collection of three files, `functions.nf`, `main.nf`, and `meta.yml`.
`functions.nf` contains utility functions for `main.nf`. The `main.nf` script
contains the process definition of the module. Finally, the `meta.yml` file
provides software metadata which can be used for input and output validation
among other things.

The first time you install an nf-core module, a `modules.json` file is also created in the root of
your workflow. This is a housekeeping file that keeps
track of your installed modules.

nf-core modules are written to a standard that makes it simpler to include them
across nf-core DSL2 workflows. However this makes it slightly more difficult to
include in a generic Nextflow DSL2 workflow.

To see how to use the nf-core modules in your workflow, follow the links below.

- [Minimal code example](./docs/MinimalCodeExample.md) - When code is easier to follow. The pages below explain in greater detail.
- [Using an nf-core module](./docs/UsingNFCoreModules.md)
- [Adapting an nf-core module](./docs/AdaptingNFCoreModules.md)
- [Making a new local module](./docs/MakingALocalModule.md)
- [Advanced module parameterisation](./docs/AdvancedParameterisation.md)
