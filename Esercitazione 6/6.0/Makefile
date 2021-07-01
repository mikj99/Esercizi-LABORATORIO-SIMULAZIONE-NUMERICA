CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ising.exe : Ising.o random.o
	$(CC) random.o Ising.o -o Ising.exe
Ising.o : Ising.cpp
	$(CC) -c Ising.cpp -o Ising.o $(CFLAGS)
random.o : random.C random.h
	$(CC) -c random.C -o random.o $(CFLAGS)
clean :
	rm *.o Ising.exe
esegui :
	./Ising.exe
