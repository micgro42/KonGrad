CC=g++
CFLAGS= -c -std=c++11



all: kongrad.o main.o
	$(CC) $(CFLAGS) kongrad.o main.o -o main

kongrad.o: kongrad.cc kongrad.hh
	$(CC) $(CFLAGS) kongrad.cc
	
main: main.o
	$(CC) main.o -o main
	
main.o: main.cc
	$(CC) $(CFLAGS) main.cc
