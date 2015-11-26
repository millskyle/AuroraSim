

files=simulation.h utility_functions.h electrostatics.cpp io.cpp electron.h

main: $(files) 
	$(CXX) -std=c++11 main.cpp -L/usr/include -lfftw3 -pthread  -orun

debug: $(files) 
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -pthread -orun && gdb ./run

valgrind: $(files) 
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -pthread -orun && valgrind ./run

clang: $(files)
	clang++ -g -std=c++11 main.cpp -L/usr/include -lfftw3 -pthread -orun

clean:
	rm run output/* energy/* rho/*

prepare:
	mkdir output energy rho
