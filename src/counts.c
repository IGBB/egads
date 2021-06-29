#include "counts.h"
#include <string.h>

/* return triangle matrix size including diagonal */
#define triangle_size(n) ((n) * ((n) + 1) / 2)
/* define triangle matrix mem alloc  */
#define trimalloc(n,t) (t*)malloc(triangle_size(n) * sizeof(t))
/* return the upper triangle index */
inline size_t get_idx(counts_t* counts, size_t a, size_t b){
    size_t i,j;
    if( a < b ){
        i = a;
        j = b;
    }else{
        i = b;
        j = a;
    }

    return (counts->n * i) + j - ((i * (i+1)) / 2);
}




counts_t* counts_init(enzyme_list_t* enzymes){
    size_t i,j;

    /* alloc mem for data struct */
    counts_t* counts = malloc(sizeof(counts_t));
    counts->n = enzymes->n;
    counts->m = triangle_size(counts->n);
    counts->d = trimalloc(enzymes->n, count_t);

    for(i = 0; i < counts->n; i++){
        for(j = i; j< counts->n; j++){
            int idx = get_idx(counts, i, j);

            /* set enzyme links */
            counts->d[idx].rare = &(enzymes->d[i]);
            counts->d[idx].freq = &(enzymes->d[j]);

            /* set all and good count array to 0  */
            memset(&(counts->d[idx].all), 0, (ALL_SIZE + 2) * sizeof(int));
            memset(&(counts->d[idx].good), 0, (GOOD_SIZE + 2) * sizeof(int));

        }
    }

    return counts;
}

counts_t* counts_site_reset(counts_t* counts){
    size_t i;
    site_t tmp = {.enz=-1, .pos=0};
    for(i = 0; i < counts->m; i++)
        counts->d[i].last = &tmp;

    return counts;
}

count_t* counts_get(counts_t* counts, size_t i, size_t j){
    size_t idx = get_idx(counts, i, j);
    return &(counts->d[idx]);
}
