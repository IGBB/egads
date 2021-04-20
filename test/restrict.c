#include "enzyme.h"
#include "sequence.h"
#include "restrict.h"
#include "data/large.test.h"
#include "data/enzymes.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


site_t output [] = {/*Enz Pos  */
                     {.enz=0,  .pos=57  },
                     {.enz=0,  .pos=314 },
                     {.enz=0,  .pos=529 },
                     {.enz=0,  .pos=655 },
                     {.enz=0,  .pos=721 },
                     {.enz=0,  .pos=763 },
                     {.enz=0,  .pos=818 },
                     {.enz=0,  .pos=913 },
                     {.enz=0,  .pos=977 },
                     {.enz=0,  .pos=1095},
                     {.enz=0,  .pos=1321},
                     {.enz=0,  .pos=1411},
                     {.enz=0,  .pos=1457},
                     {.enz=0,  .pos=1714},
                     {.enz=0,  .pos=2235},
                     {.enz=0,  .pos=2335},
                     {.enz=0,  .pos=2398}
};

int main(){
    size_t length;
    enzyme_list_t enzymes = {.n=enzyme_list_length, .m=3, .d=enzyme};

        length = strlen(large_test_fa);


        site_list_t * sites = site_list_init();

        sites = restrict_scan(&enzymes, large_test_fa, length, sites);

        int i;
        printf("  Start     End  Strand Enzyme_name Restriction_site  5prime"
               "  3prime  5frag  3frag 5primerev 3primerev 5fragrev 3fragrev\n");

        for(i=0; i<sites->n; i++){
            uint64_t bases =  enzymes.d[sites->d[i].enz].pattern;
            char rest_site [16];
            int j;
            for(j = 0; j < 16; j++){
                rest_site[15-j] = seq_char[bases&0xF];
                bases >>= 4;
            }
            rest_site[enzymes.d[sites->d[i].enz].length] = 0;
            printf("%7d %7d %7c %-11s %-16s %7d %7d"
                   "      .      .         .         .        .        .\n",
                   sites->d[i].pos+1,
                   sites->d[i].pos + enzymes.d[sites->d[i].enz].length ,
                   '+', enzymes.d[sites->d[i].enz].name,
                   rest_site,
                   sites->d[i].pos + enzymes.d[sites->d[i].enz].c[0],
                   sites->d[i].pos + enzymes.d[sites->d[i].enz].c[1]);
        }
    return 0;
}
