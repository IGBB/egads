#include "args.h"

#include "klib/ketopt.h"
#include <string.h>
#include <stdlib.h>

const char* const help_message =
  "Usage: egads [OPTION...] --rare <rare[,...]> --freq <freq[,...]>\n"
  "                         --genome <fasta[.gz]> > <output.html>\n"
  "egads -- Electronically Guided Digestion Selection \n\n"
  "  -r, --rare <rare[,...]>    Comma delimited list of rare cutting enzymes\n"
  "  -f, --freq <freq[,...]>    Comma delimited list of frequent cutting\n"
  "                               enzymes (can overlap with list of rare for\n"
  "                               single enzyme digestion results)\n"
  "  -g, --genome <fasta[.gz]>  FASTA file of sequences to digest (can be gzip\n"
  "                               compressed)\n"
  "  -h, --html FILE            HTML output file (default: stdout)\n"
  "  -b, --bed FILE             BED output file (tab-delimited file with \n"
  "                               columns: sequnece id, start, end, \n"
  "                               rare-freq enzyme names, length of fragment)\n"
  "  -t, --title title          Title of html report (default: genome file name)\n"
  "  -?, --help                 Give this help list\n"
  "Report bugs to github.com/IGBB/egads.\n";


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
    
    { "help", ko_no_argument, '?' },
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
  while ((c = ketopt(&opt, argc, argv, 1, "r:f:g:h:b:t:?", longopts)) >= 0) {
    switch(c){
      case 'g': arguments.genome  = opt.arg; break;
      case 't': arguments.title  = opt.arg; break;

      case 'r':
        arguments.rare.n=parse_enzymes(opt.arg, arguments.rare.d);
        if(arguments.rare.n == 0) {
          fprintf(stderr, "Cannot parse rare enzyme list");
          fprintf(stderr, help_message);
          exit(EXIT_FAILURE);
        }
        break;
      case 'f':
        arguments.freq.n=parse_enzymes(opt.arg, arguments.freq.d);
        if(arguments.freq.n == 0) {
          fprintf(stderr, "Cannot parse freq enzyme list");
          fprintf(stderr, help_message);
          exit(EXIT_FAILURE);
        }
        break;
      case 'b':
      case 'h':
        tmp = fopen(opt.arg, "w");

        if(c == 'b') arguments.bed = tmp;
        if(c == 'h') arguments.html = tmp;

        if(tmp == NULL){
          perror("Cannot open output file\n");
          fprintf(stderr, help_message);
          exit(EXIT_FAILURE);
        }

        break;

      case '?':
       printf(help_message);
       exit(EXIT_SUCCESS);
    };
  }

  /* If no title is given, use genome filename */
  if(arguments.title == NULL)
    arguments.title = arguments.genome;

  /* If no output arguments are given, then default to svg on stdout */
  if(arguments.html == NULL && arguments.bed == NULL)
    arguments.html = stdout;


  if(arguments.genome == NULL ||
     arguments.rare.n == 0 ||
     arguments.freq.n == 0) {
          fprintf(stderr, "The arguments --genome, --rare, and --freq "
                  "must be given\n");
          fprintf(stderr, help_message);
          exit(EXIT_FAILURE);
  }

  return arguments;
}
