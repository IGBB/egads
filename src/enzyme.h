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
        buffer_string[300],/* \0 delimited list of buffers */
        supplier[20];      /* supplier list */

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

/*******************************************************************************
 * enzyme_is_compat
 *
 *   1. enzyme_t*         :: enzyme 1
 *   2. enzyme_t*         :: enzyme 2
 *   r. int               :: boolean
 *
 *   return 1 if the same string in the buffer list is found, else 0. A special
 *   case was added to return 2 if both enzymes are supplied by NEB but don't
 *   have the same buffer name.
 *
 ******************************************************************************/

int enzyme_is_compat(enzyme_t*, enzyme_t*);


#endif // __ENZYME_H_
