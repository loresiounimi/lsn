CFLAGS = -Wall -O3 -std=c++14

all : r.o random.o statistics.o
	g++ random.o r.o statistics.o -o r

r.o : r.cpp
	g++ -c r.cpp -o r.o $(CFLAGS)

random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)

statistics.o : statistics.cpp statistics.h
	g++ -c statistics.cpp -o statistics.o $(CFLAGS)