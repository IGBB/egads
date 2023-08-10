#include <stdlib.h>
#include <stdio.h>
#include "enzyme.h"
#include "counts.h"

#define INCBIN_PREFIX incbin_
#define INCBIN_STYLE INCBIN_STYLE_SNAKE
#define INCBIN_SILENCE_BITCODE_WARNING
#include "incbin/incbin.h"


/* html/d3.v5.min.js */
INCBIN(d3, "html/d3.v5.min.js");

/* html/draw.d3.css */
INCBIN(css, "html/draw.d3.css");

/* html/draw.table.js */
INCBIN(table, "html/draw.table.js");

/* html/draw.ridge.js */
INCBIN(ridge, "html/draw.ridge.js");

/* html/draw.gel.js */
INCBIN(gel, "html/draw.gel.js");

/* html/draw.html */
INCBIN(html, "html/draw.html");

void print_html (FILE* out,
                 char* name,
                 counts_t* counts,
                 long genome_size,
                 int mutation_frags){

    size_t i, j ;


    fprintf( out,
             "<!doctype html>\n"
             "<html>\n"
             "    <head>\n"
             "        <meta charset=\"utf-8\">\n"
             "        <title>%s</title>\n"
             "        <meta name=\"description\" content=\"\">\n"
             "        <meta name=\"viewport\" \n"
             "              content=\"width=device-width, initial-scale=1\">\n",
             name);

    /* print html/draw.d3.css */
    fprintf(out, "<style>\n");
    fwrite(incbin_css_data, sizeof(char), incbin_css_size, out);
    fprintf(out, "</style>\n");
    /* print html/draw.table.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(incbin_d3_data, sizeof(char), incbin_d3_size, out);
    fprintf(out, "</script>\n");
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(incbin_table_data, sizeof(char), incbin_table_size, out);
    fprintf(out, "</script>\n");
    /* print html/draw.ridge.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(incbin_ridge_data, sizeof(char), incbin_ridge_size, out);
    fprintf(out, "</script>\n");
    /* print html/draw.gel.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(incbin_gel_data, sizeof(char), incbin_gel_size, out);
    fprintf(out, "</script>\n");

    /* json output */
    fprintf(out, "<script type=\"application/json\" id=\"data\">\n");
    fprintf(out, "[");
    for(i = 0; i < counts->m; i++){
        if(i > 0) fprintf(out, ",");

        fprintf(out, "{");

        fprintf(out, "\"name\": \"%s - %s\",",
               counts->d[i].rare->name,
               counts->d[i].freq->name);

        /* warn if either enzymes are blunt-ended */
        fprintf(out, "\"blunt\": %d,",
                counts->d[i].rare->blunt | counts->d[i].freq->blunt);

        /* warn if enzyme pairs don't have same buffer */
        fprintf(out, "\"compat\": %d,",
                enzyme_is_compat(counts->d[i].rare, counts->d[i].freq));


        fprintf(out, "\"genome_size\": %ld,", genome_size);
        fprintf(out, "\"mutation\": %d,", mutation_frags);

        fprintf(out, "\"all\": [");
        for(j = 0; j < ALL_SIZE+2; j++){
            if(j > 0) fprintf(out, ",");
            fprintf(out, "%d", counts->d[i].all[j]);
        }
        fprintf(out, "],");

        fprintf(out, "\"good\": [");
        for(j = 0; j < GOOD_SIZE+2; j++){
            if(j > 0) fprintf(out, ",");
            fprintf(out, "%d", counts->d[i].good[j]);
        }

        fprintf(out, "]");
        fprintf(out, "}");
    }

    fprintf(out, "]");

    fprintf(out, "\n</script>\n");
    fprintf(out, "</head>\n");
    fwrite(incbin_html_data, sizeof(char), incbin_html_size, out);
    fprintf( out, "</html>" );
}
