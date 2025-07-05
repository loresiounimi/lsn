CFLAGS = -Wall -O3 --std=c++14

esercizio2.2.1 : esercizio2.2.1.o random.o 
	g++ random.o esercizio2.2.1.o  -o esercizio2.2.1
esercizio2.2.2 : esercizio2.2.2.o random.o 
	g++ random.o esercizio2.2.2.o -o esercizio2.2.2
esercizio2.2.1.o : esercizio2.2.1.cpp 
	g++ -c esercizio2.2.1.cpp -o esercizio2.2.1.o $(CFLAGS)
esercizio2.2.2.o : esercizio2.2.2.cpp 
	g++ -c esercizio2.2.2.cpp -o esercizio2.2.2.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
