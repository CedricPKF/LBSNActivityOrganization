CFLAGS = -g -Wall
CFLAGS += -O3

compute_SMK: compute_SMK.o functions.o
	g++ -pthread -o $@ compute_SMK.o functions.o $(CFLAGS)

compute_SMK.o: compute_SMK.cpp functions.h
	g++ -c compute_SMK.cpp $(CFLAGS)

functions.o: functions.cpp functions.h Node.h
	g++ -c functions.cpp $(CFLAGS)

clean:
	rm -rf compute_SMK *.o
