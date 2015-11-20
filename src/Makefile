

main: data_structure.h project_viewpoint.cpp utility_functions.h
	$(CXX) -std=c++11 main.cpp -L/usr/include -lfftw3 -orun
debug: data_structure.h project_viewpoint.cpp utility_functions.h
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clang: data_structure.h project_viewpoint.cpp utility_functions.h
	clang++ -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clean:
	rm run
