CXXFLAGS = -Wall --pedantic -g -std=c++11
LIBS = `root-config --libs`
INCS = `root-config --cflags`
PARAL = -fopenmp


main.x: main.C random.o
	g++ -o main.x main.C random.o ${CXXFLAGS} ${LIBS} ${INCS} ${PARAL}

random.o : random.C random.h
	g++ -c random.C -o random.o ${CXXFLAGS} ${LIBS} ${INCS} ${PARAL}

compila: main.x

esegui:
	./main.x

clean:
	rm *.o
	rm *.x

c++11: 
	scl enable devtoolset-8 bash




