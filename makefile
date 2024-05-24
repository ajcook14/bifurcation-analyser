all: main

main.o: main.cpp main.h solvers.h neuron.h capd_error.h
	g++ -pg -g main.cpp -c -o main.o -O2 `capd-config --cflags`

small_funcs.o: small_funcs.cpp small_funcs.h
	g++ -pg -g small_funcs.cpp -c -o small_funcs.o -O2 `capd-config --cflags`

solvers.o: solvers.cpp solvers.h small_funcs.h capd_error.h
	g++ -pg -g solvers.cpp -c -o solvers.o -O2 `capd-config --cflags`

neuron.o: neuron.cpp neuron.h
	g++ -pg -g neuron.cpp -c -o neuron.o -O2 `capd-config --cflags`

main: main.o small_funcs.o solvers.o neuron.o
	g++ -pg -g $? -o $@ `capd-config --libs`

test: capd_test.cpp
	g++ $? -o $@ -O2 `capd-config --cflags --libs`

clean:
	rm -f *.o main test

