CC=g++
CFLAGS= -std=c++11 -Wall -pedantic -g -O0 -DBOOST_ALL_DYN_LINK -I/users/stud/micgro42/boost/include/
LDFLAGS= 
LDLIBS= -L/users/stud/micgro42/boost/lib/ -lpthread -lboost_log -lboost_unit_test_framework


all: kongrad.o main.o unittest
	$(CC) $(LDFLAGS) kongrad.o main.o $(LDLIBS) -o main

kongrad.o: kongrad.cc kongrad.hh
	$(CC) $(CFLAGS) -c kongrad.cc
	
main.o: main.cc
	$(CC) $(CFLAGS) -c main.cc
	
unittest: test.o kongrad.o
	$(CC) $(LDFLAGS) test.o kongrad.o $(LDLIBS) -o unittest

test.o: test.cc
	$(CC) $(CFLAGS) -c test.cc

clean:
	rm -f *.o main unittest
