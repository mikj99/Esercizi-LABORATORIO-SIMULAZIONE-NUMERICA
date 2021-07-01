PARAL = -fopenmp


main.x: main.C 
	g++ -o main.x main.C ${PARAL}

random.o : MolDyn_NVE.C MolDyn_NVE.h
	g++ -c MolDyn_NVE.C -o MolDyn_NVE.o ${CXXFLAGS} ${LIBS} ${INCS} ${PARAL}

compila: main.x

esegui:
	./main.x

clean:
	rm *.o
	rm *.x

c++11: 
	scl enable devtoolset-8 bash




