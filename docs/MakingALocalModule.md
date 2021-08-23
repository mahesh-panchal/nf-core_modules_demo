## Making A Local Module

When a module is not available from nf-core, it's time to
make your own locally. The best way to start is to use
```bash
nf-core modules create <tool/subtool>
```
which then instructs the `nf-core` tool to make a new module
template in `modules/local/<tool>_<subtool>.nf`.
