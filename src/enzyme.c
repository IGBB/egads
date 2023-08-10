#include "enzyme.h"
#include "sequence.h"

#include <string.h>
#include <stdlib.h>

#define INCBIN_PREFIX incbin_
#define INCBIN_STYLE INCBIN_STYLE_SNAKE
#define INCBIN_SILENCE_BITCODE_WARNING
#include "incbin/incbin.h"

INCTXT(msbuffmin, "msbuffmin.txt");

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


enzyme_list_t* load_enzymes(char* input, char** names, int len){
    int i;
    char *s = NULL;
    size_t slen = 0;
    ssize_t rlen;

    /* Open enzyme file. If file isn't provided, open mem stream to saved
     * restriction enzymes */
    FILE* file;
    if(input != NULL){
        file = fopen(input, "r");
        /* die if file can't be opened */
        if( file == NULL ){
            perror("Can't open enzyme file");
            exit(1);
        }
    }else{
        file = fmemopen((void *) incbin_msbuffmin_data, incbin_msbuffmin_size, "r");
    }

    /* pattern string to be encoded. NOTE: Must be equal to the number of bases
     * the enzyme_t pattern element can hold */
    char pattern[16];

    enzyme_list_t * list = enzyme_list_init();

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
        }

        /* Add space for new enzyme */
        e = enzyme_list_push(list);

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

        /* Suppliers */
        strncpy(e->supplier, strtok2(NULL,';'), 20);


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

    fclose(file);

    return list;
};


/* rudementary buffer compatibility check. Doesn't work well since the only
 * buffer data from REBASE is incomplete.  */
int enzyme_is_compat(enzyme_t* a, enzyme_t* b){
    int i = 0, j = 0;

    /* loop until end of a or b buffer list is reached  */
    while(i < a->buffer_length && j < b->buffer_length){

        /* compare current buffer names */
        int compare = strcmp(a->buffer_list[i], b->buffer_list[j]);

        /* if the two buffers are the same, return true. Else, increment the
         * index of the smaller value ( works since the lists are sorted ) */
        if(compare == 0 ) return 1;
        else if ( compare < 0 ) i++;
        else if ( compare > 0 ) j++;
    }

    /* Since NEB is moving many of their enzymes to the smartcut buffer but not
     * updating the msbufmin file, a special case is added to return "maybe" if
     * both enzymes are supplied by NEB */
    if(strchr(a->supplier, 'N') != NULL &&
       strchr(b->supplier, 'N') != NULL){
        return 2;
    }


    return 0;
}
