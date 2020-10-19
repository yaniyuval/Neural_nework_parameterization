/*
 system-function for Fortran
 (C) Marat Khairoutdinov */

#include <stdio.h>

#ifdef __cplusplus 
extern "C" {
#endif

void systemf(const char *string) {system(string);}

void systemf_(const char *string) {system(string);}


#ifdef __cplusplus
}
#endif
