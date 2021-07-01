CXXFLAGS = -Wall --pedantic -g -std=c++11

main.x: main.C random.o genetic.o
	mpicxx -o main.x main.C random.o genetic.o ${CXXFLAGS}

random.o : random.C random.h
	g++ -c random.C -o random.o ${CXXFLAGS}

genetic.o : genetic.C genetic.h 
	g++ -c genetic.C -o genetic.o ${CXXFLAGS}


compila: main.x

esegui:
	mpirun -np 4 main.x

clean:
	rm *.o
	rm *.x

c++11: 
	scl enable devtoolset-8 bash

loadmpi:
	module load mpi


