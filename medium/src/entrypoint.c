// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_rowSum_extendr(void *dll);

void R_init_rowSum(void *dll) {
    R_init_rowSum_extendr(dll);
}
