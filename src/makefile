##
# ddrad viz
#
# @file
# @version 0.1

CC = gcc
CFLAGS = -Wall -O3  -lm -ggdb
PROG=../egads

.PHONY: all clean

all: $(PROG)

clean:
	rm -f *.o $(PROG)

html/d3.v5.min.js:
	wget -O $@ https://d3js.org/d3.v5.min.js
html/matter.min.css:
	wget -O $@ https://github.com/finnhvman/matter/releases/download/0.2.2/matter.min.css

msbuffmin.txt:
	wget -O $@ ftp://ftp.neb.com/pub/rebase/msbuffmin.txt


enzyme.o : sequence.o msbuffmin.txt
restrict.o : enzyme.o sequence.o
html.o : html/d3.v5.min.js html/draw.d3.css html/draw.gel.js html/draw.html html/draw.ridge.js html/draw.table.js html/matter.min.css
main.o : restrict.o enzyme.o sequence.o html.o


$(PROG): main.o restrict.o enzyme.o sequence.o counts.o html.o args.o klib/kstring.o
	$(CC) -o $@ $^ $(CFLAGS) -lz
# end
