CC=g++
CFLAGS= -c -std=c++11



all: kongrad.cc kongrad.hh main.cc
	gcc -o $(CFLAGS) main kongrad.cc main.cc
	
main: main.o
	gcc main.o -o main
	
main.o: main.cc
	$(CC) $(CFLAGS) main.cc