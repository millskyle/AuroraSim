

files=data_structure.h project_viewpoint.cpp utility_functions.h electrostatics.cpp io.cpp electron.h

main: $(files) 
	$(CXX) -std=c++11 main.cpp -L/usr/include -lfftw3 -orun
debug: data_structure.h project_viewpoint.cpp utility_functions.h electrostatics.cpp
	$(CXX) -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clang: data_structure.h project_viewpoint.cpp utility_functions.h electrostatics.cpp
	clang++ -g -std=c++11 main.cpp -L/usr/include -lfftw3 -orun

clean:
	rm run
