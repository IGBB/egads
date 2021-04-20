#include "enzyme.h"
#include "sequence.h"

#include <string.h>
#include <stdlib.h>

enzyme_list_t* enzyme_list_init(){
    enzyme_list_t * r = malloc(sizeof(enzyme_list_t));
    r->m = 0;
    r->n = 0;
    r->d = NULL;
    return r;
}

inline enzyme_t* enzyme_list_push(enzyme_list_t* list){
    enzyme_t * ret;

    if(list->n >= list->m){
        list->m += 100;
        list->d = realloc(list->d, list->m*sizeof(enzyme_t));
    }

    ret = &(list->d[list->n]);
    list->n++;

    return ret;

}

enzyme_list_t* load_enzymes(FILE* input, char** names, int len){
    int i;
    char *s = NULL;
    size_t slen = 0;
    ssize_t rlen;

    /* if input isn't provided, open mem stream to saved restriction enzymes */
    FILE* file = input;
    if(file == NULL)
        file = fmemopen(_binary_emboss_e_txt_start,
                        _binary_emboss_e_txt_size,
                        "r");


    /* pattern string to be encoded. NOTE: Must be equal to the number of bases
     * the enzyme_t pattern element can hold */
    char pattern[16];

    enzyme_list_t * list = enzyme_list_init();

    /* if names are given, adjust the size of the enzyme list to match the
     * length of the given names; so we can add the enzymes in the correct
     * order */
    if(len){
        list->n = len-1;
        enzyme_list_push(list);
    }

    /* read each line in enzyme database */
    while((rlen = getline(&s, &slen, file)) != -1 ){

        /* skip comments */
        if(s[0] == '#') continue;

        /* get name (first tab column) */
        char * name = strtok(s, "\t");

        enzyme_t * e;

        /* If not reading all, check if this enzyme matches any fo the names */
        if(len > 0){
            /* loop through names, stopping when found a match */
            for(i = 0;
                i < len && strcmp(names[i], name) != 0 ;
                i++);
            /* if i == len then no match found; skip entry */
            if(i == len) continue;

            /* set entry point to correct location */
            e = &(list->d[i]);

        /* else enzyme is added to end of list */
        }else {
            e = enzyme_list_push(list);
        }


       strncpy(e->name, name, 32);

       /* clear pattern and copy pattern string to array*/
       memset(pattern , 0 , 16);
       strncpy(pattern, strtok(NULL,"\t"), 16);

       e->length = atoi(strtok(NULL,"\t"));
       e->ncuts  = atoi(strtok(NULL,"\t"));
       e->blunt  = atoi(strtok(NULL,"\t"));

       for(i = 0; i < 4; i++)
           e->c[i]  = atoi(strtok(NULL,"\t"));

       /* Encode pattern   */
       uint64_t base;
       e->pattern = 0;  /* set pattern to 0, use ~or~ (|) to store each base */
       for(i = 0; i < 16; i ++){
           base = (uint64_t) pattern[i]; /* extract base from pattern */
           base = seq_table[base] ;      /* encode base to 4bit flag */
           base <<= ((15-i) * 4);        /* shift base to correct position - 1st
                                          * base is the most significat 4 bits*/
           e->pattern |= base;           /* store encoded base in pattern */
       }

    }


    /* if file is mem buffer, close file */
    if(input == NULL) fclose(file);

    return list;
};
