#include "args.h"

#include "klib/ketopt.h"
#include <string.h>
#include <stdlib.h>

size_t parse_enzymes(char * string, char * list []){
    size_t len = 0;

    list[len] = strtok(string, ",");
    while(list[len] != NULL){
        len++;
        list[len] = strtok(NULL, ",");
    }

    return len;
}



static ko_longopt_t longopts[] = {
    { "rare",  ko_required_argument, 'r' },
    { "freq",  ko_required_argument, 'f' },

    { "genome", ko_required_argument, 'g' },
    { "title", ko_required_argument, 't' },

    { "html", ko_required_argument, 'h' },
    { "bed", ko_required_argument, 'b' },

    {NULL, 0, 0}
  };


arguments_t parse_options(int argc, char **argv) {
  arguments_t arguments = {
                                .rare.n = 0,
                                .freq.n = 0,

                                .genome = NULL,

                                .html  = NULL,
                                .title = NULL,
                                .bed   = NULL
  };


  ketopt_t opt = KETOPT_INIT;

  int  c;
  FILE* tmp;
  while ((c = ketopt(&opt, argc, argv, 1, "r:f:g:h:b:t", longopts)) >= 0) {
    switch(c){
      case 'g': arguments.genome  = opt.arg; break;
      case 't': arguments.title  = opt.arg; break;

      case 'r': arguments.rare.n=parse_enzymes(opt.arg, arguments.rare.d); break;
      case 'f': arguments.freq.n=parse_enzymes(opt.arg, arguments.freq.d); break;

      case 'b':
      case 'h':
        tmp = fopen(opt.arg, "w");

        if(c == 'b') arguments.bed = tmp;
        if(c == 'h') arguments.html = tmp;

        if(tmp == NULL){
          perror("Cannot open output file\n");
          /* perror(help_message);   */
          exit(EXIT_FAILURE);
        }

        break;

    };
  }

  /* If no title is given, use genome filename */
  if(arguments.title == NULL)
    arguments.title = arguments.genome;

  /* If no output arguments are given, then default to svg on stdout */
  if(arguments.html == NULL && arguments.bed == NULL)
    arguments.html = stdout;

    return arguments;
}
