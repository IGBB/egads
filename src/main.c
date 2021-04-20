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


extern void print_html (FILE*, char*, enzyme_list_t*, counts_t*, size_t, size_t, long, int);


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
    FILE* emboss_e = NULL;

    if(file != NULL){
        emboss_e = fopen(file, "r");

        if(emboss_e == NULL){
            perror("Can't open emboss_e.txt file");
            exit(1);
        }
    }
    enzyme_list_t * enzymes = load_enzymes(emboss_e, names, len);

    if(file != NULL)
        fclose(emboss_e);

    return enzymes;
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

    for(i =0; i < enzymes->n; i++){
        fprintf(stderr, "%s, ", enzymes->d[i].name);
    }
    fprintf(stderr, "\n");



    site_t ** matrix =  malloc( (rare_length * freq_length) * sizeof(site_t*) );
    counts_t * counts =  malloc( (rare_length * freq_length) * sizeof(counts_t) );

    memset(counts, 0, (rare_length * freq_length) * sizeof(counts_t));

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
        for(i = 0; i < rare_length * freq_length; i++){
            matrix[i] = &tmp;
        }


        /* loop through the sites, collating them to rare-freq pairs  */
        for(i = 0; i < sites->n; i++){
            site_t* cur = &(sites->d[i]);

            size_t start, end, stride;
            /* TODO: Since we expect freq cutter to be seen more than the rare,
             * it may be benificial to make rare cutters columns and freq cutter
             * rows*/
            /* rare cutter (row access)*/
            if(cur->enz < rare_length){
                start  = cur->enz * freq_length;
                end    = start + freq_length;
                stride = 1;
            /* freq cutter (column access)*/
            } else {
                start  = cur->enz - rare_length;
                end    = rare_length*freq_length;
                stride = freq_length;
            }


            /* compare new and old sites */
            for(j = start; j < end ; j+=stride){

                /* get length of current fragment, adjust if longer than 10kb */
                int len = cur->pos - matrix[j]->pos;
                len = (len > ALL_SIZE)? ALL_SIZE+1 : len;



                counts[j].all[len]++;

                /* adjust length again if longer than 1kb */
                len = (len > GOOD_SIZE)? GOOD_SIZE+1 : len;
                /* store fragment if between different  */
                if( matrix[j]->enz != cur->enz &&
                                matrix[j]->enz != -1 )
                    counts[j].good[len]++;

                /* replace previous site */
                matrix[j]=cur;
            }


        }

        /* final fragments */
        for(i = 0; i < rare_length * freq_length; i++){

            /* get length of last fragment, adjust if longer than 10kb */
                int len =  seq->seq.l - matrix[i]->pos;
                len = (len > ALL_SIZE)? ALL_SIZE+1 : len;

                counts[i].all[len]++;

        }



        /* clear the sites vector  */
        site_list_clear(sites);

    }
    kseq_destroy(seq);
    gzclose(fp);

    /* free the sites vector */
    site_list_free(sites);
    sites = NULL;

    print_html(stdout,
               argv[3],
               enzymes,
               counts,
               rare_length,
               freq_length,
               genome_size,
               mutation_frags);


    return 0;
}
