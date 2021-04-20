#include "enzyme.h"

#include <stdlib.h>
#include <stdio.h>

#include "data/enzymes.h"


int main(int argc, char** argv){

    FILE* emboss_e = fopen("../data/emboss_e.txt", "r");
    if(emboss_e == NULL){
        perror("Can't open emboss_e.txt file");
        return 1;
    }

    enzyme_list_t * list = load_enzymes(emboss_e,
                                        enzyme_names,
                                        enzyme_list_length);


    if(list->n != enzyme_list_length){
        printf("Failed to reads correct number of enzymes: FAIL\n");
        return 1;
    }

    int i ;
    for(i = 0; i < list->n; i++){
        if(strcmp(list->d[i].name, enzyme[i].name) != 0){
            printf("Enzymes are not in correct order: \n"
                   "    Expected: %s   Observed: %s\n",
                   enzyme[i].name, list->d[i].name);
            return 1;
        }

        if(list->d[i].pattern != enzyme[i].pattern ||
           list->d[i].length  != enzyme[i].length  ||
           list->d[i].ncuts   != enzyme[i].ncuts   ||
           list->d[i].blunt   != enzyme[i].blunt   ){
            printf("Enzyme %s does not have expected value:\n");
            printf("    %s\t%llx\t%d\t%d\t%d\n",
                   list->d[i].name,
                   list->d[i].pattern,
                   list->d[i].length,
                   list->d[i].ncuts,
                   list->d[i].blunt);
            printf("    %s\t%llx\t%d\t%d\t%d\n",
                   enzyme[i].name,
                   enzyme[i].pattern,
                   enzyme[i].length,
                   enzyme[i].ncuts,
                   enzyme[i].blunt);
            return 1;
        }
    }

    return 0;

}
