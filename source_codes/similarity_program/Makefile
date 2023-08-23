CFLAGS = -g -Wall
CFLAGS += -O3

alg: alg.o functions.o
	g++ -pthread -o $@ alg.o functions.o $(CFLAGS)

alg.o: alg.cpp functions.h
	g++ -c alg.cpp $(CFLAGS)

functions.o: functions.cpp functions.h Node.h
	g++ -c functions.cpp $(CFLAGS)

clean:
	rm -rf alg *.o
