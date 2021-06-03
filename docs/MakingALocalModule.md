## Making A Local Module

When a module is not available from nf-core, it's time to
make your own locally. The best way to start is to use
```bash
nf-core modules create -a <github_author> -t <software/tool> .
```

This requires the presence of a `main.nf` in the current directory,
which then instructs the `nf-core` tool to make a new module
template in `modules/local/<software>/<tool>.nf`.
