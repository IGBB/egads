#ifndef ARGS_H_
#define ARGS_H_
#include <stdio.h>

extern const char* const program_version;

typedef struct {
    struct { size_t n; char * d[64]; }
        rare, freq;
    char * genome, * title;
    FILE * html, * bed;
} arguments_t;

arguments_t parse_options(int argc, char **argv);


#endif // ARGS_H_
