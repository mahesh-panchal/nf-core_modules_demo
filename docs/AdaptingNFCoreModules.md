## Adapting an nf-core module

Not all modules can be used as they are. While in some
cases modifying the `args` parameter in the `module['software']`
Map is enough, some need to
be modified for use in a workflow. A common example is the
MultiQC module in which the input declaration is often
customised to the workflow.

nf-core DSL2 workflows store adapted modules to the path
`modules/local`. A copy of the `functions.nf` file needs to be copied
here, and the `main.nf` has been renamed to something
based on the tool name. e.g. `multiqc_illumina.nf`.

Let's follow this convention too. Begin by installing the module as before.
```bash
nf-core modules install multiqc
```
which can be found in `modules/nf-core/modules/multiqc/`.
Then make the folder `modules/local`, and copy the `functions.nf` script
there.
```bash
mkdir -p modules/local
cp modules/nf-core/modules/multiqc/functions.nf modules/local/functions.nf
```
Next make a copy of the `main.nf` renaming the script after the tool.
```bash
cp modules/nf-core/modules/multiqc/main.nf modules/local/multiqc.nf
```

If you want, you can then remove the installed module using
```bash
nf-core modules remove multiqc
```

The next step is to change the necessary parts of the code
to fit your workflow, for example the input declaration.

Let's make a custom module for MultiQC
```bash
nf-core modules install multiqc
mkdir -p modules/local
cp modules/nf-core/modules/multiqc/functions.nf modules/local/functions.nf
cp modules/nf-core/modules/multiqc/main.nf modules/local/multiqc.nf
```
and then change the input from:
```nextflow
    input:
    path multiqc_files
```
to:
```nextflow
    path multiqc_config
    path software_versions
    path ('fastqc/*')
```

We can then include the module into our workflow:
```nextflow
include { MULTIQC } from './modules/local/multiqc' addParams(options:modules['multiqc'])
```
and add :
```
params {
    modules {
        'fastqc' {
            <fastqc module params>
        }
        'multiqc' {
            publish_dir = 'MultiQC_Report'
        }
    }
}
```
to your `conf/modules.config` file.
