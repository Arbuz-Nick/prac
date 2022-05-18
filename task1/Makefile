all: main

main: *.cpp *.h
	g++ *.cpp -o prog -std=c++17 -fopenmp

omp_polus: *.cpp *.h
	xlc++ *.cpp -o prog -Wall -std=c++11 -qsmp=omp -fopenmp

clean:
	rm -rf ./prog
