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


/* mimic strtok function; but don't skip empty delimeters */
char* strtok2(char* string, char token){
    /* setup locally scoped but persistant variable */
    static char* copy;

    /* set persistant varible if string is set (new string to tokenize)  */
    if(string != NULL) copy = string;

    /* if finished tokenizing string, return null */
    if(copy == NULL) return NULL;

    char* ret = copy;

    /* find token in string */
    char* tok = strchr(copy, token);

    /* if token is found */
    if(tok != NULL){
        /* replace token with '\0' */
        *tok = '\0';
        /* set persistant variable to the remaining string*/
        copy = tok + 1;
    } else {
        /*clear persistant variable*/
        copy = NULL;
    }

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
        file = fmemopen(_binary_msbuffmin_txt_start,
                        _binary_msbuffmin_txt_size,
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

        /* skip lines that begin with space */
        if(s[0] == ' ') continue;

        /* get name (first column) */
        char * name = strtok2(s, ';');

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


        /* Safely copy name to the enzyme struct  */
        strncpy(e->name, name, 32);

        /* clear pattern and copy pattern string to array*/
        memset(pattern , 0 , 16);
        strncpy(pattern, strtok2(NULL,';'), 16);
        e->length = strlen(pattern);

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

        /* 5` and 3` cut site  */
        e->c[0] = atoi(strtok2(NULL,';'));
        e->c[1] = atoi(strtok2(NULL,';'));

        /* Enzyme is blunt if cut sites are equal */
        e->blunt = (e->c[0] == e->c[1]);

        /* Methylation data (discarded) */
        strtok2(NULL,';');

        /* Enzyme Type (discarded) */
        strtok2(NULL,';');

        /* Suppliers (discarded) */
        strtok2(NULL,';');


        /* Buffers (semicolon terminated, comma delimited list) NOTE: some of
         * the buffer names contain semicolons, e.g "O:Toyobo NheI; RsaI; SfiI"
         */

        /* Safely copy the rest of the line to buffer_string */
        strncpy(e->buffer_string, strtok2(NULL,'\n'), 300);
        /* Ensure safe string */
        e->buffer_string[299] = '\0';
        /* Remove semicolon termination */
        e->buffer_string[strlen(e->buffer_string) - 1 ] = '\0';

        /* Tokenize buffer_string, sorting and storing locations in
         * buffer_list */
        char* key = strtok2(e->buffer_string, ',');
        e->buffer_length = 0;

        /* Loop until no more buffers found */
        while(key != NULL){

            /* find correct spot int buffer_list (insertion sort) */
            for(i = e->buffer_length-1;
                i >= 0 && strcmp(key, e->buffer_list[i]) == 1;
                i--)
                e->buffer_list[i+1] = e->buffer_list[i];

            /* insert key, and increment length */
            e->buffer_list[i+1] = key;
            e->buffer_length += 1;

            /* read next buffer */
            key = strtok2(NULL, ',');
        }

    }

    free(s);

    /* if file is mem buffer, close file */
    if(input == NULL) fclose(file);

    return list;
};
