#ifndef __ENZYME_H_
#define __ENZYME_H_

#include <stdint.h>
#include <stdio.h>

/* Saved data from msbuffmin.txt at compile time */
extern const char _binary_msbuffmin_txt_start[];
extern const char _binary_msbuffmin_txt_end[];
#define _binary_msbuffmin_txt_size                              \
    (_binary_msbuffmin_txt_end - _binary_msbuffmin_txt_start)




typedef struct {
    char name [32],        /* Name of Enzyme */
        *buffer_list[20],  /* Sorted list of buffer names */
        buffer_string[300];/* \0 delimeted list of buffers */

    uint64_t pattern;      /* Encodded Recognition site, able to hold 16 bases,
                            * needs to be right padded with 0*/
    int length,            /* length of enzyme */
        blunt,             /* true if blunt end cut, false if sticky */
        c[2],              /* 5` and 3` cut site */
        buffer_length;     /*  */
} enzyme_t;

typedef struct {
    size_t m,n;
    enzyme_t *d;
} enzyme_list_t;

/*******************************************************************************
 * load_enzymes
 *
 *   1. FILE*          :: enzyme file, NULL cause the file to be read from object
 *                           cached at compile time
 *   2. char**         :: list of enzyme names to load
 *   3. int            :: number of name in list, 0 if all are loaded
 *   r. enzyme_list_t* :: returned enzyme list
 *
 *   Load enzymes from msbuffmin.txt file.
 ******************************************************************************/
enzyme_list_t* load_enzymes(FILE*, char**, int );

#endif // __ENZYME_H_
