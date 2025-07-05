CC = g++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo -lstdc++fs

all: simulator_forward.exe simulator_backward.exe
simulator_forward.exe : system.o NSL_SIMULATOR_forward.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR_forward.o -o simulator_forward.exe $(AFLAGS)
NSL_SIMULATOR_forward.o : NSL_SIMULATOR_forward.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR_forward.cpp -o NSL_SIMULATOR_forward.o $(CFLAGS)
simulator_backward.exe : system.o NSL_SIMULATOR_backward.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR_backward.o -o simulator_backward.exe $(AFLAGS)
NSL_SIMULATOR_backward.o : NSL_SIMULATOR_backward.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR_backward.cpp -o NSL_SIMULATOR_backward.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
