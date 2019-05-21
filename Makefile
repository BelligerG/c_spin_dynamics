CXXFLAGS = -Wall -Werror -I/home/links/cms242/eigen-eigen-b3f3d4950030 -std=c++11 -fopenmp -Ofast

spinDynamics: main.o spin_dynamics_plusplus.o
	g++ -fopenmp -o test main.o spin_dynamics_plusplus.o
	./test

