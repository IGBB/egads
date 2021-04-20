#include "sequence.h"

#include <string.h>
#include <ctype.h>
#include <stdio.h>

const char * sequence[] = { "garkbdctymvhu", "ctymvhgarkbda", "ccccccccccga" };
const char * rev_comp[] = { "adbkraghvmytc", "thvmytcdbkrag", "tcgggggggggg" };

int main (){
    int test, i, len;
    char tmp [128];

    /* test encode and decode */
    printf("Validate Encode -> Decode \n");
    printf("  -  Test 0 should FAIL becuase of the u -> t encoding error\n");
    for(test = 0; test < 3; test++){
        len = strlen(sequence[test]);
        for(i = 0; i < len; i++)
            tmp[i] = tolower(seq_char[seq_table[sequence[test][i]]]);
        tmp[len] = '\0';

        printf("\tTest %d: %16s == %-16s : ", test, sequence[test], tmp);
        if(strcmp(sequence[test], tmp) == 0)
            printf("PASS\n");
        else
            printf("FAIL\n");

    }

    printf("Validate reverse complement encode \n");
    for(test = 0; test < 3; test++){
        len = strlen(sequence[test]);
        for(i = 0; i < len; i++)
            tmp[len - i - 1] = tolower(seq_char[comp_table[sequence[test][i]]]);
        tmp[len] = '\0';

        printf("\tTest %d: %16s == %-16s : ", test, rev_comp[test], tmp);
        if(strcmp(rev_comp[test], tmp) == 0)
            printf("PASS\n");
        else
            printf("FAIL\n");

    }


    return 0;
}
