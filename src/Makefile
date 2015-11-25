

files=simulation.h utility_functions.h electrostatics.cpp io.cpp electron.h

main: $(files) 
	$(CXX) -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

debug: $(files) 
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun && gdb ./run

valgrind: $(files) 
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun && valgrind ./run

clang: $(files)
	clang++ -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clean:
	rm run
