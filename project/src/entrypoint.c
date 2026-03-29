#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_rust_extendr(void *dll);

void R_init_rust(DllInfo *dll) {
    R_init_rust_extendr(dll);
}
