CFLAGS = -Wall -O3 --std=c++11

all : main.o random.o
	g++ random.o main.o -o main
main.o : main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)