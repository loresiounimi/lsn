CFLAGS = -Wall -O3 --std=c++11

all : chiquadro.o random.o statistics.o
	g++ random.o chiquadro.o statistics.o -o chiquadro
chiquadro.o : chiquadro.cpp
	g++ -c chiquadro.cpp -o chiquadro.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
statistics.o : statistics.cpp statistics.h
	g++ -c statistics.cpp -o statistics.o $(CFLAGS)
