

main: data_structure.h project_viewpoint.cpp utility_functions.h
	$(CXX) -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clang: data_structure.h project_viewpoint.cpp utility_functions.h
	clang -std=c++11 main.cpp -orun

clean:
	rm run
