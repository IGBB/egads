#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <math.h>

#include "enzyme.h"
#include "restrict.h"
#include "klib/kseq.h"
#include "counts.h"



KSEQ_INIT(gzFile, gzread);


extern void print_html (FILE*, char*, counts_t*, size_t, long, int);


int parse_enzymes(char * string, char * list []){
    int len = 0;

    list[len] = strtok(string, ",");
    while(list[len] != NULL){
        len++;
        list[len] = strtok(NULL, ",");
    }

    return len;
}

enzyme_list_t * get_enzymes(char * file, char ** names, int len){
    FILE* msbuffmin = NULL;

    if(file != NULL){
        msbuffmin = fopen(file, "r");

        if(msbuffmin == NULL){
            perror("Can't open msbuffmin.txt file");
            exit(1);
        }
    }
    enzyme_list_t * enzymes = load_enzymes(msbuffmin, names, len);

    if(file != NULL)
        fclose(msbuffmin);

    return enzymes;
}

/* return triangle matrix size including diagonal */
#define triangle_size(n) ((n) * ((n) + 1) / 2)
/* return the upper triangle index */
inline size_t get_idx(size_t n, size_t a, size_t b){
    size_t i,j;
    if( a < b ){
        i = a;
        j = b;
    }else{
        i = b;
        j = a;
    }

    return (n * i) + j - ((i * (i+1)) / 2);
}

int main(int argc, char **argv) {
    size_t i,j;
    char * enzyme_file = NULL;
    gzFile fp;
    kseq_t *seq;
    int l;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <rare[,...]> <freq[,...]> <in.fasta>\n", argv[0]);
        return 1;
    }


    char *names[128];
    /* load rare cutters into name list */
    size_t rare_length = parse_enzymes(argv[1], names);
    /* load freq cutters into name list */
    size_t freq_length = parse_enzymes(argv[2], &(names[rare_length]));

    /* get enzyme definitions from file */
    enzyme_list_t * enzymes = get_enzymes(enzyme_file,
                                          names,
                                          rare_length + freq_length);


    // Create triangle matrix to store last site seen for each enzyme pair
    site_t ** matrix = malloc(triangle_size(enzymes->n) * sizeof(counts_t) );

    // Create triangle matrix to store counts of each enzyme pair, setting all to 0
    counts_t * counts =  malloc(triangle_size(enzymes->n) * sizeof(counts_t) );
    memset(counts, 0, triangle_size(enzymes->n) * sizeof(counts_t) );
    for(i = 0; i < enzymes->n; i++){
        for(j = i; j< enzymes->n; j++){
            int idx = get_idx(enzymes->n, i, j);
            counts[idx].rare = &(enzymes->d[i]);
            counts[idx].freq = &(enzymes->d[j]);
       }
    }

    site_list_t * sites = site_list_init();

    long genome_size = 0;
    int mutation_frags = 0;

    fp = gzopen(argv[3], "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {

        /* Calculate the mutation rate for current sequence
         *        The regression for cellular organisms is
         *        −0.81 + 0.68log10(G), with r2 = 0.80.
         *     https://doi.org/10.1016/j.tig.2010.05.003
         * */
        mutation_frags += round(-0.81 + 0.68 * log10((double) seq->seq.l));


        /* Calculate genome size */
        genome_size += seq->seq.l;

        /* find all sites that are cut by the given enzymes */
        sites = restrict_scan(enzymes, seq->seq.s, seq->seq.l, sites);

        /* clear matrix for new sequence */
        site_t tmp = {.enz=-1, .pos=0};
        for(i = 0; i < triangle_size(enzymes->n); i++)
            matrix[i] = &tmp;

        /* loop through the sites, collating them to rare-freq pairs  */
        for(i = 0; i < sites->n; i++){
            site_t* cur = &(sites->d[i]);

            /* compare new and old sites */
            for(j = 0; j < enzymes->n; j++){
                int idx = get_idx(enzymes->n, cur->enz, j);

                /* get length of current fragment, adjust if longer than 10kb */
                int len = cur->pos - matrix[idx]->pos;
                len = (len > ALL_SIZE)? ALL_SIZE+1 : len;

                counts[idx].all[len]++;

                /* adjust length again if longer than 1kb */
                len = (len > GOOD_SIZE)? GOOD_SIZE+1 : len;

                /* store fragment if between different enzymes, or if current
                 * pair is single enzyme digestion */
                if((matrix[idx]->enz != -1        &&
                    matrix[idx]->enz != cur->enz) ||
                    j == cur->enz)
                    counts[idx].good[len]++;

                /* replace previous site */
                matrix[idx]=cur;
            }
        }

        /* final fragments */
        for(i = 0; i < enzymes->n; i++){
            for(j = 0; j < enzymes->n; j++){
                int idx = get_idx(enzymes->n,i,j);
            /* get length of last fragment, adjust if longer than 10kb */
                int len =  seq->seq.l - matrix[idx]->pos;
                len = (len > ALL_SIZE)? ALL_SIZE+1 : len;

                counts[idx].all[len]++;
            }
        }

        /* clear the sites vector  */
        site_list_clear(sites);

    }
    kseq_destroy(seq);
    gzclose(fp);

    /* free the sites vector */
    site_list_free(sites);
    sites = NULL;

    // free data structures no longer needed
    free(matrix);

    size_t count_size = 0;
    for(i = 0; i < triangle_size(enzymes->n); i++){
        int rare = 0, freq = 0;
        // Check rare list for the two enzymes
        for(j = 0; j < rare_length; j++){
            if(strcmp(counts[i].rare->name, names[j]) == 0)
                rare |= 1;
            if(strcmp(counts[i].freq->name, names[j]) == 0)
                freq |= 1;
        }

        // Check freq list for the two enzymes
        for(j = rare_length; j < rare_length+freq_length; j++){
            if(strcmp(counts[i].rare->name, names[j]) == 0)
                rare |= 2;
            if(strcmp(counts[i].freq->name, names[j]) == 0)
                freq |= 2;
        }

        /* skip unless both enzymes in the ueser given list */
        if(freq == 0 || rare == 0 || (freq | rare) != 3 ) continue;


        /* swap enzymes if rare is in freq list or vise-versa */
        if(rare == 2 || freq == 1){
            enzyme_t* swp = counts[i].rare;
            counts[i].rare = counts[i].freq;
            counts[i].freq = swp;
        }

        /* move current counts to next available space */
        counts[count_size] = counts[i];
        count_size++;

    }


    print_html(stdout,
               argv[3],
               counts,
               count_size,
               genome_size,
               mutation_frags);


    return 0;
}
