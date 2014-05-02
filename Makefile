CC=g++
CFLAGS= -std=c++11 -Wall -pedantic -g -O0 -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/
LDFLAGS= 
LDLIBS= -L/users/stud/micgro42/boost/lib -lpthread -lboost_log -lboost_unit_test_framework -lboost_thread -lboost_filesystem -lboost_date_time -lboost_chrono -lboost_system
SRCDIR=src/

all: kongrad.o main.o sparseMatrix.o unittest
	$(CC) $(LDFLAGS) kongrad.o sparseMatrix.o main.o $(LDLIBS) -o main

kongrad.o: $(SRCDIR)kongrad.cc $(SRCDIR)kongrad.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)kongrad.cc
	
sparseMatrix.o: $(SRCDIR)sparseMatrix.cc $(SRCDIR)sparseMatrix.hh
	$(CC) $(CFLAGS) -c $(SRCDIR)sparseMatrix.cc
	
main.o: $(SRCDIR)main.cc
	$(CC) $(CFLAGS) -c $(SRCDIR)main.cc
	
unittest: test.o kongrad.o sparseMatrix.o
	$(CC) $(LDFLAGS) test.o kongrad.o sparseMatrix.o $(LDLIBS) -o unittest

test.o: $(SRCDIR)test.cc
	$(CC) $(CFLAGS) -c $(SRCDIR)test.cc

clean:
	rm -f *.o main unittest
