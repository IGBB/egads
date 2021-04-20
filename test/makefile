##
# unit test
#
# @file
# @version 0.1
CC = gcc
CFLAGS = -Wall -O3 -I../src/ -ftree-vectorize -fopt-info-vec -fopt-info-vec-missed
TESTS=sequence restrict enzyme

.PHONY: all clean $(TESTS)

all: $(TESTS)

$(TESTS) : % : %.test
	./$<

enzyme.test : enzyme.c ../src/sequence.o ../src/enzyme.o
sequence.test : sequence.c ../src/sequence.o
restrict.test : restrict.c ../src/sequence.o ../src/restrict.o

%.test :
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^

../src/%.o :
	 $(MAKE) $(MFLAGS) -C ../src/ $(notdir $@)


clean:
	rm -f *.o



# end
