CFLAGS = -Wall -O3 --std=c++11

all : main.o random.o statistics.o
	g++ random.o main.o statistics.o -o main
main.o : main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
statistics.o : statistics.cpp statistics.h
	g++ -c statistics.cpp -o statistics.o $(CFLAGS)