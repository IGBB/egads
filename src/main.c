#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <math.h>

#include "enzyme.h"
#include "restrict.h"
#include "klib/kseq.h"
#include "counts.h"
#include "args.h"


KSEQ_INIT(gzFile, gzread);


extern void print_html (FILE*, char*, counts_t*, long, int);


int main(int argc, char **argv) {
    size_t i,j;
    char * enzyme_file = NULL;
    gzFile fp;
    kseq_t *seq;
    int l;

    arguments_t arguments = parse_options(argc, argv);

    struct {size_t n; char * d[128]; } names = {.n=0};
    for(i = 0; i < arguments.rare.n; i++){
        names.d[names.n] = arguments.rare.d[i];
        names.n++;
    }
    for(i = 0; i < arguments.freq.n; i++){
        names.d[names.n] = arguments.freq.d[i];
        names.n++;
    }

    /* get enzyme definitions from file */
    enzyme_list_t * enzymes = load_enzymes(enzyme_file,
                                          names.d, names.n);

    counts_t * counts = counts_init(enzymes);
    site_list_t * sites = site_list_init();

    long genome_size = 0;
    int mutation_frags = 0;

    fp = gzopen(arguments.genome, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {


        /* Calculate genome size */
        genome_size += seq->seq.l;

        /* find all sites that are cut by the given enzymes */
        sites = restrict_scan(enzymes, seq->seq.s, seq->seq.l, sites);

        /* clear last sites for new sequence */
        counts_site_reset(counts);

        /* loop through the sites, collating them to rare-freq pairs  */
        for(i = 0; i < sites->n; i++){
            site_t* cur = &(sites->d[i]);

            /* compare new and old sites */
            for(j = 0; j < enzymes->n; j++){
                count_t * c = counts_get(counts, cur->enz, j);

                /* Skip all checks for first site on contig */
                if(c->last != NULL){

                    /* Print bed formatted hit to stderr
                     *TODO: make a cmd flag to toggle this */
                    if(arguments.bed != NULL && c->last->enz != cur->enz)
                        fprintf(arguments.bed, "%s\t%u\t%u\t%s-%s\t%u\n",
                                seq->name.s,  c->last->pos, cur->pos+1,
                                c->rare->name, c->freq->name,
                                cur->pos - c->last->pos);


                    /* get length of current fragment, adjust if longer than 10kb */
                    int len = cur->pos - c->last->pos;
                    if (len > ALL_SIZE) len = ALL_SIZE+1;

                    c->all[len]++;

                    /* adjust length again if longer than 1kb */
                    if (len > GOOD_SIZE) len = GOOD_SIZE+1;

                    /* store fragment if between different enzymes, or if current
                     * pair is single enzyme digestion */
                    if(c->last->enz != cur->enz ||
                       (ssize_t)j == cur->enz)
                        c->good[len]++;
                }



                /* replace previous site */
                c->last=cur;
            }
        }

        /* final fragments */
        for(i = 0; i < counts->m; i++){

            /* Skip section if neither enzyme was found on current contig */
            if(counts->d[i].last == NULL) continue;

            /* get length of last fragment, adjust if longer than 10kb */
            int len =  seq->seq.l - counts->d[i].last->pos;
            len = (len > ALL_SIZE)? ALL_SIZE+1 : len;

            counts->d[i].all[len]++;
        }

        /* clear the sites vector  */
        site_list_clear(sites);

    }
    kseq_destroy(seq);
    gzclose(fp);

    /* free the sites vector */
    site_list_free(sites);
    sites = NULL;

    size_t count_size = 0;
    for(i = 0; i < counts->m; i++){
        int rare = 0, freq = 0;
        // Check rare list for the two enzymes
        for(j = 0; j < arguments.rare.n; j++){
            if(strcmp(counts->d[i].rare->name, arguments.rare.d[j]) == 0)
                rare |= 1;
            if(strcmp(counts->d[i].freq->name, arguments.rare.d[j]) == 0)
                freq |= 1;
        }

        // Check freq list for the two enzymes
        for(j = 0; j < arguments.freq.n; j++){
            if(strcmp(counts->d[i].rare->name, arguments.freq.d[j]) == 0)
                rare |= 2;
            if(strcmp(counts->d[i].freq->name, arguments.freq.d[j]) == 0)
                freq |= 2;
        }

        /* skip unless both enzymes in the ueser given list */
        if(freq == 0 || rare == 0 || (freq | rare) != 3 ) continue;


        /* swap enzymes if rare is in freq list or vise-versa */
        if(rare == 2 || freq == 1){
            enzyme_t* swp = counts->d[i].rare;
            counts->d[i].rare = counts->d[i].freq;
            counts->d[i].freq = swp;
        }

        /* move current counts to next available space */
        counts->d[count_size] = counts->d[i];
        count_size++;

    }
    counts->m = count_size;

    /* Calculate the mutation rate for current sequence
     *        The regression for cellular organisms is
     *        −0.81 + 0.68log10(G), with r2 = 0.80.
     *     https://doi.org/10.1016/j.tig.2010.05.003
     * */
    mutation_frags = round(-0.81 + 0.68 *
                           log10((double) genome_size ));
    fprintf(stderr, "genome_size = %ld\nmutation rate = %lf\n",
            genome_size,
            -0.81 + 0.68 * log10((double) genome_size/1000000.0 ));


    if(arguments.html != NULL)
        print_html(arguments.html,
                   arguments.title,
                   counts,
                   genome_size,
                   mutation_frags);


    return 0;
}
