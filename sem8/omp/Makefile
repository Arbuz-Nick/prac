all: main

main_no_barriers: *.cpp *.h
	g++ main_no_barriers.cpp matrix.cpp -o prog -std=c++11 -fopenmp

main_barriers: *.cpp *.h
	g++ main_barriers.cpp matrix.cpp -o prog -std=c++11 -fopenmp

omp_polus_no_barriers: *.cpp *.h
	xlc++ main_no_barriers.cpp matrix.cpp -o prog -Wall -std=c++11 -qsmp=omp -fopenmp

omp_polus_barriers: *.cpp *.h
	xlc++ main_barriers.cpp matrix.cpp -o prog -Wall -std=c++11 -qsmp=omp -fopenmp

clean:
	rm -rf ./prog
