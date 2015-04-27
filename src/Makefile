all: tmm

util.o: util.c util.h
	g++ -c -g -Wall util.c

poleCompute.o: poleCompute.c util.h
	g++ -c -g -Wall poleCompute.c

residueCompute.o: residueCompute.c util.h
	g++ -c -g -Wall residueCompute.c

tmm.o: tmm.c util.h
	g++ -c -g -Wall tmm.c

main.o: main.c util.h
	g++ -c -g -Wall main.c

tmm: util.o poleCompute.o residueCompute.o tmm.o main.o
	g++ util.o poleCompute.o residueCompute.o tmm.o main.o -o tmm

clean: 
	rm -rf *.o tmm *~ *.core
