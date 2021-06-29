#ifndef __COUNTS_H_
#define __COUNTS_H_

#include <stdlib.h>
#include "enzyme.h"
#include "restrict.h"

#define ALL_SIZE  10000
#define GOOD_SIZE  1000



typedef struct {
    int all  [ALL_SIZE +2 ],
        good [GOOD_SIZE + 2];
    enzyme_t *rare, *freq;
    site_t* last;
} count_t;

typedef struct {
    size_t n,m;
    count_t* d;
} counts_t;

// Create triangle matrix to store counts of each enzyme pair
counts_t* counts_init(enzyme_list_t*);
// Set last cut site to 0
counts_t* counts_site_reset(counts_t*);
// Get count data for enzyme pair
count_t* counts_get(counts_t*, size_t, size_t);
#endif // __COUNTS_H_
