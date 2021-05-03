#include "restrict.h"

#include <stdlib.h>

#include "sequence.h"

site_list_t* site_list_init(){
    site_list_t * r = malloc(sizeof(site_list_t));
    r->n = 0;

    /* set an initial size so the push logic is simpler */
    r->m = 100;
    r->d = malloc(r->m * sizeof(site_t));

    return r;
}

/* clear list without free'ing */
void site_list_clear(site_list_t* list){
    list->n=0;
}

/* free list data pointer and list pointer */
void site_list_free(site_list_t* list){
    free(list->d);
    free(list);
}


inline site_t* site_list_push(site_list_t* list){
    site_t * ret;

    if(list->n == list->m){
        list->m *= 2;
        list->d = realloc(list->d, list->m*sizeof(site_t));
    }

    ret = &(list->d[list->n]);
    list->n++;

    return ret;

}



inline void restrict_kmer(enzyme_list_t* enzymes,
                             uint64_t kmer,
                             int pos,
                             site_list_t* list){
    size_t i,j;
    for(i = 0; i < enzymes->n; i++){
        uint64_t test = kmer & enzymes->d[i].pattern;
        int acc = 0;

        /* loop through each base of the 16-mer and count the matches */
        for(j = 0; j < 16; j++){
            /* test if j-th base matched pattern */
            acc += ((test & (0xFllu << (j*4))) != 0);
        }

        /* if enzyme matches current kmer */
        if(acc == enzymes->d[i].length) {
            site_t* cur = site_list_push(list);
            cur->enz = i;
            cur->pos = pos;
        }
    }
}


site_list_t*  restrict_scan(enzyme_list_t* enzymes, char* seq, int len, site_list_t* list){
    int i;
    uint64_t kmer = 0;

    /* fill 16-mer */
    for(i = 0; i < 15; i++)
        kmer = (kmer << 4) + seq_table_strict[(int) seq[i]];

    /* check 16-mers against enzyme pattern */
    for(; i < len; i++){
        kmer = (kmer << 4) + seq_table_strict[(int) seq[i]];
        restrict_kmer(enzymes, kmer, i-15, list);
    }

    /* empty 16-mer */
    for(i = 0; i < 15; i++){
        kmer <<= 4;
        restrict_kmer(enzymes, kmer, len-15+i, list);
    }


    return list;
};
