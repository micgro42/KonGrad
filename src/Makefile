CC=g++
CFLAGS= -std=c++11 -Wall -pedantic -g -O0



all: kongrad.o main.o test
	$(CC) kongrad.o main.o -o main

kongrad.o: kongrad.cc kongrad.hh
	$(CC) $(CFLAGS) -c kongrad.cc
	
main.o: main.cc
	$(CC) $(CFLAGS) -c main.cc
	
test: test.o
	$(CC) test.o -o test

test.o: test.cc
	$(CC) $(CFLAGS) -c test.cc

clean:
	rm -f *.o main