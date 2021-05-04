#ifndef __COUNTS_H_
#define __COUNTS_H_

#include "enzyme.h"

#define ALL_SIZE  10000
#define GOOD_SIZE  1000


typedef struct {
    int all  [ALL_SIZE +2 ],
        good [GOOD_SIZE + 2];
    enzyme_t *rare, *freq;
} counts_t;

#endif // __COUNTS_H_
