#g++ -O2 bisection.cpp -o bisection `capd-config --cflags --libs`
#-frounding-math -D__USE_FILIB__ -frounding-math -I/usr/local/include -L/usr/local/lib -lcapd -lboost_filesystem -lboost_system -lboost_regex
all: main

main.o: main.cpp solvers.h neuron.h
	g++ -O2 main.cpp -c -o main.o -frounding-math -D__USE_FILIB__ -I/usr/local/include

small_funcs.o: small_funcs.cpp small_funcs.h
	g++ -O2 small_funcs.cpp -c -o small_funcs.o -frounding-math -D__USE_FILIB__ -I/usr/local/include

solvers.o: solvers.cpp solvers.h small_funcs.h
	g++ -O2 solvers.cpp -c -o solvers.o -frounding-math -D__USE_FILIB__ -I/usr/local/include

neuron.o: neuron.cpp neuron.h
	g++ -O2 neuron.cpp -c -o neuron.o -frounding-math -D__USE_FILIB__ -I/usr/local/include

main: main.o small_funcs.o solvers.o neuron.o
	g++ $? -o $@ -L/usr/local/lib -lcapd -lboost_filesystem -lboost_system -lboost_regex

test: capd_test.cpp
	g++ -O2 $? -o $@ -frounding-math -D__USE_FILIB__ -frounding-math -I/usr/local/include -L/usr/local/lib -lcapd -lboost_filesystem -lboost_system -lboost_regex

clean:
	rm -f *.o main test

