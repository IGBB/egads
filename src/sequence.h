#ifndef __SEQUENCE_H_
#define __SEQUENCE_H_

/* encode to char */
extern const char seq_char[16];

/* char -> encoding tables form ACTG and ambiguous bases, everything else is 0  */
extern const unsigned char seq_table[256];



/* complement(char) -> encoding tables form ACTG and ambiguous bases, everything
 * else is 0 */
extern const unsigned char comp_table[256];

/* Same as seq_table but N don't match anything  */
extern const unsigned char seq_table_strict[256];

#endif // __SEQUENCE_H_
