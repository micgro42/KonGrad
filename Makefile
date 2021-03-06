### Compile settings
CC=g++
CXXWARNINGS=-Wall -Wextra -pedantic -Wshadow -Wlogical-op
CFLAGS= -Ofast -fopenmp -std=c++11 $(CXXWARNINGS) -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/

### Linker settings
LDFLAGS= -fopenmp
LDLIBS= -L/users/stud/micgro42/boost/lib -lpthread -lboost_log -lboost_unit_test_framework -lboost_thread -lboost_filesystem -lboost_date_time -lboost_chrono -lboost_system

### Debug Settings
DBGDIR = Debug
DBGCFLAGS = -g3 -DDEBUG

### Release Settings
RLSDIR = Release
RELCFLAGS = -DNDEBUG

SRCDIR=src/

all: unittest main

debug: CFLAGS += $(DBGCFLAGS)
debug: OUTDIR = $(DBGDIR)
debug: all

main: linsysequ.o main.o geom_pbc.o timedif.o eivtris.o
	$(CC) $(LDFLAGS) eivtris.o timedif.o linsysequ.o main.o $(LDLIBS) -o main

linsysequ.o: $(SRCDIR)linsysequ.cc $(SRCDIR)linsysequ.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)linsysequ.cc
	
sparseMatrix.o: $(SRCDIR)sparseMatrix.cc $(SRCDIR)sparseMatrix.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)sparseMatrix.cc
	
main.o: $(SRCDIR)main.cc 
	$(CC) $(CFLAGS) -c $(SRCDIR)main.cc
	
unittest: test.o linsysequ.o geom_pbc.o timedif.o eivtris.o
	$(CC) $(LDFLAGS) eivtris.o timedif.o linsysequ.o test.o $(LDLIBS) -o unittest

test.o: $(SRCDIR)test.cc linsysequ.o
	$(CC) $(CFLAGS) -c $(SRCDIR)test.cc

geom_pbc.o: $(SRCDIR)geom_pbc.c
	$(CC) $(CFLAGS) -c $(SRCDIR)geom_pbc.c

timedif.o: $(SRCDIR)timedif.h $(SRCDIR)timedif.c
	$(CC) $(CFLAGS) -c $(SRCDIR)timedif.c
	
eivtris.o: $(SRCDIR)eivtris.h $(SRCDIR)eivtris.c
	$(CC) $(CFLAGS) -fpermissive -c $(SRCDIR)eivtris.c

clean:
	rm -f *.o main unittest
	
.PHONY: all clean debug release 
