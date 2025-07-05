CFLAGS = -Wall -O3 --std=c++14

all : main.o random.o integral.o
	g++ random.o main.o integral.o -o main
main.o : main.cpp funzioni.h
	g++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
integral.o : integral.cpp integral.h funzioni.h
	g++ -c integral.cpp -o integral.o $(CFLAGS)