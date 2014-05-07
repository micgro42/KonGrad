CC=g++
#CFLAGS= -fopenmp -std=c++11 -Wall -pedantic -g -O0 -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/ 
CFLAGS= -fopenmp -std=c++11 -O3 -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/
#CFLAGS= -std=c++11 -O3 -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/
LDFLAGS= -fopenmp
LDLIBS= -L/users/stud/micgro42/boost/lib -lpthread -lboost_log -lboost_unit_test_framework -lboost_thread -lboost_filesystem -lboost_date_time -lboost_chrono -lboost_system
SRCDIR=src/

all: unittest main

main: kongrad.o main.o geom_pbc.o timedif.o
	$(CC) $(LDFLAGS) timedif.o kongrad.o main.o $(LDLIBS) -o main

kongrad.o: $(SRCDIR)kongrad.cc $(SRCDIR)kongrad.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)kongrad.cc
	
sparseMatrix.o: $(SRCDIR)sparseMatrix.cc $(SRCDIR)sparseMatrix.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)sparseMatrix.cc
	
main.o: $(SRCDIR)main.cc 
	$(CC) $(CFLAGS) -c $(SRCDIR)main.cc
	
unittest: test.o kongrad.o geom_pbc.o timedif.o
	$(CC) $(LDFLAGS) timedif.o kongrad.o test.o $(LDLIBS) -o unittest

test.o: $(SRCDIR)test.cc kongrad.o
	$(CC) $(CFLAGS) -c $(SRCDIR)test.cc

geom_pbc.o: $(SRCDIR)geom_pbc.c
	$(CC) $(CFLAGS) -c $(SRCDIR)geom_pbc.c

timedif.o: $(SRCDIR)timedif.h $(SRCDIR)timedif.c
	$(CC) $(CFLAGS) -c $(SRCDIR)timedif.c

clean:
	rm -f *.o main unittest
	

