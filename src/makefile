##
# ddrad viz
#
# @file
# @version 0.1

CC = gcc
CFLAGS = -Wall -O3 -lm -ggdb
PROG=../egads

.PHONY: all clean

all: $(PROG)

clean:
	rm -f *.o $(PROG)

html/d3.v5.min.js:
	wget -O $@ https://d3js.org/d3.v5.min.js

emboss_e.txt:
	wget -O $@ ftp://ftp.neb.com/pub/rebase/emboss_e.txt


d3.bin.o:    html/d3.v5.min.js
css.bin.o:   html/draw.d3.css
html.bin.o:  html/draw.html
gel.bin.o:   html/draw.gel.js
ridge.bin.o: html/draw.ridge.js
table.bin.o: html/draw.table.js

emboss.bin.o : emboss_e.txt

%.bin.o:
	ld -b binary -r -o $@ $<

html.o : d3.bin.o css.bin.o html.bin.o gel.bin.o ridge.bin.o table.bin.o
enzyme.o : sequence.o emboss.bin.o
restrict.o : enzyme.o sequence.o
main.o : restrict.o enzyme.o sequence.o html.o


$(PROG): main.o restrict.o enzyme.o sequence.o html.o d3.bin.o css.bin.o \
			html.bin.o gel.bin.o ridge.bin.o table.bin.o emboss.bin.o
	$(CC) -o $@ $(CFLAGS) -lz $^
# end
