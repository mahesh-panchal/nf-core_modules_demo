## Adapting an nf-core module

Not all modules can be used as is. Some need to
be modified for use in a workflow. A common example is the
MultiQC module in which the input declaration is often
customised to the workflow.

nf-core DSL2 workflows store adapted modules to the path
`modules/local`. A copy of the `functions.nf` file resides
here, and the `main.nf` has been renamed to something
based on the tool name. e.g. `multiqc_illumina.nf`.

Let's follow this convention too. Begin by installing the module as before.
```bash
nf-core modules install --tool multiqc .
```
which can be found in `modules/nf-core/software/multiqc/`.
Then make the folder `modules/local`, and copy the `functions.nf` script
there.
```bash
mkdir -p modules/local
cp modules/nf-core/software/multiqc/functions.nf modules/local/functions.nf
```
Next make a copy of the `main.nf` renaming the script after the tool.
```bash
cp modules/nf-core/software/multiqc/main.nf modules/local/multiqc.nf
```

If you want, you can then remove the installed module using
```bash
nf-core modules remove --tool multiqc .
```

The next step is to change the necessary parts of the code
to fit your workflow, for example the input declaration.
