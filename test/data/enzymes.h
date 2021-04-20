#ifndef __ENZYMES_H_
#define __ENZYMES_H_
#include "enzyme.h"

int enzyme_list_length = 3;
char *enzyme_names[] = { "AluI", "HindII", "HindIII"};

enzyme_t enzyme[] = {
                       { /* AluI  AGCT  4  2  1   2 2 0 0 */
                       .name = "AluI",
                       .pattern=0x1428000000000000,
                       .length=4,
                       .ncuts=2,
                       .blunt=1,
                       .c={2,2,0,0}
                       },
                       { /*  HindII  GTYRAC  6  2  1  3 3 0 0*/
                       .name = "HindII",
                       .pattern=0x48A5120000000000,
                       .length=6,
                       .ncuts=2,
                       .blunt=1,
                       .c={3,3,0,0}
                       },
                       { /* HindIII  AAGCTT  6  2  0  1 5 0 0 */
                       .name = "HindIII",
                       .pattern=0x1142880000000000,
                       .length=6,
                       .ncuts=2,
                       .blunt=0,
                       .c={1,5,0,0}
                       }
};


#endif // __ENZYMES_H_
