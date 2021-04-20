#ifndef __RESTRICT_H_
#define __RESTRICT_H_

#include "enzyme.h"

typedef struct  {
    ssize_t enz;
    unsigned int pos;
} site_t;

typedef struct {
    size_t n,m;
    site_t *d;
} site_list_t;


site_list_t* site_list_init();
void site_list_clear(site_list_t*);
void site_list_free(site_list_t*);

/* scans the sequence for interactions with the given enzymes. Returns a sorted
 * vector of restriction sites */
site_list_t* restrict_scan(enzyme_list_t*, char*, int, site_list_t* );

#endif // __RESTRICT_H_
