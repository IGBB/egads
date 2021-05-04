#include <stdlib.h>
#include <stdio.h>
#include "enzyme.h"
#include "counts.h"

/* html/d3.v5.min.js */
extern const char _binary_html_d3_v5_min_js_start[];
extern const char _binary_html_d3_v5_min_js_end[];
#define _binary_html_d3_v5_min_js_size                              \
    (_binary_html_d3_v5_min_js_end - _binary_html_d3_v5_min_js_start)

/* html/draw.d3.css */
extern const char _binary_html_draw_d3_css_start[];
extern const char _binary_html_draw_d3_css_end[];
#define _binary_html_draw_d3_css_size                              \
    (_binary_html_draw_d3_css_end - _binary_html_draw_d3_css_start)

/* html/draw.table.js */
extern const char _binary_html_draw_table_js_start[];
extern const char _binary_html_draw_table_js_end[];
#define _binary_html_draw_table_js_size                              \
    (_binary_html_draw_table_js_end - _binary_html_draw_table_js_start)

/* html/draw.ridge.js */
extern const char _binary_html_draw_ridge_js_start[];
extern const char _binary_html_draw_ridge_js_end[];
#define _binary_html_draw_ridge_js_size                              \
    (_binary_html_draw_ridge_js_end - _binary_html_draw_ridge_js_start)

/* html/draw.gel.js */
extern const char _binary_html_draw_gel_js_start[];
extern const char _binary_html_draw_gel_js_end[];
#define _binary_html_draw_gel_js_size                              \
    (_binary_html_draw_gel_js_end - _binary_html_draw_gel_js_start)

/* html/draw.html */
extern const char _binary_html_draw_html_start[];
extern const char _binary_html_draw_html_end[];
#define _binary_html_draw_html_size                              \
    (_binary_html_draw_html_end - _binary_html_draw_html_start)



void print_html (FILE* out,
                 char* name,
                 counts_t* counts,
                 size_t length,
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
    fwrite(_binary_html_draw_d3_css_start, sizeof(char),
           _binary_html_draw_d3_css_size, out);
    fprintf(out, "</style>\n");
    /* print html/draw.table.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(_binary_html_d3_v5_min_js_start, sizeof(char),
           _binary_html_d3_v5_min_js_size, out);
    fprintf(out, "</script>\n");
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(_binary_html_draw_table_js_start, sizeof(char),
           _binary_html_draw_table_js_size, out);
    fprintf(out, "</script>\n");
    /* print html/draw.ridge.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(_binary_html_draw_ridge_js_start, sizeof(char),
           _binary_html_draw_ridge_js_size, out);
    fprintf(out, "</script>\n");
    /* print html/draw.gel.js */
    fprintf(out, "<script type=\"text/javascript\">\n");
    fwrite(_binary_html_draw_gel_js_start, sizeof(char),
           _binary_html_draw_gel_js_size, out);
    fprintf(out, "</script>\n");

    /* json output */
    fprintf(out, "<script type=\"application/json\" id=\"data\">\n");
    fprintf(out, "[");
    for(i = 0; i < length; i++){
        if(i > 0) fprintf(out, ",");

        fprintf(out, "{");

        fprintf(out, "\"name\": \"%s - %s\",",
               counts[i].rare->name,
               counts[i].freq->name);

        /* warn if either enzymes are blunt-ended */
        fprintf(out, "\"blunt\": %d,",
                counts[i].rare->blunt | counts[i].freq->blunt);

        /* warn if enzyme pairs don't have same buffer */
        fprintf(out, "\"compat\": %d,",
                enzyme_is_compat(counts[i].rare, counts[i].freq));


        fprintf(out, "\"genome_size\": %ld,", genome_size);
        fprintf(out, "\"mutation\": %d,", mutation_frags);

        fprintf(out, "\"all\": [");
        for(j = 0; j < ALL_SIZE+2; j++){
            if(j > 0) fprintf(out, ",");
            fprintf(out, "%d", counts[i].all[j]);
        }
        fprintf(out, "],");

        fprintf(out, "\"good\": [");
        for(j = 0; j < GOOD_SIZE+2; j++){
            if(j > 0) fprintf(out, ",");
            fprintf(out, "%d", counts[i].good[j]);
        }

        fprintf(out, "]");
        fprintf(out, "}");
    }

    fprintf(out, "]");

    fprintf(out, "\n</script>\n");
    fprintf(out, "</head>\n");


    fwrite(_binary_html_draw_html_start, sizeof(char),
           _binary_html_draw_html_size, out);



    fprintf( out, "</html>" );
}
