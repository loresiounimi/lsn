CFLAGS = -Wall -O3 --std=c++14
LDFLAGS = -larmadillo

all : main.o random.o chromosomes.o individuo.o
	g++ random.o main.o chromosomes.o individuo.o -o main $(LDFLAGS)
main.o : main.cpp 
	g++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
chromosomes.o : chromosomes.cpp chromosomes.h
	g++ -c chromosomes.cpp -o chromosomes.o $(CFLAGS)
individuo.o : individuo.cpp individuo.h
	g++ -c individuo.cpp -o individuo.o $(CFLAGS)

