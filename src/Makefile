CC=g++
CFLAGS= -c -std=c++11 -Wall -pedantic -g -O0



all: kongrad.o main.o
	$(CC) kongrad.o main.o -o main

kongrad.o: kongrad.cc kongrad.hh
	$(CC) $(CFLAGS) kongrad.cc
	
main.o: main.cc
	$(CC) $(CFLAGS) main.cc

clean:
	rm -f *.o main