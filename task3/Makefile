all: omp_polus
#	bsub < ./lsf_files/100/4.lsf
main: *.cpp *.h
	g++ *.cpp -o prog -std=c++17 -fopenmp

omp_polus: *.cpp *.h
	xlc++ main.cpp CG_copy.cpp matrix.cpp -o prog -Wall -std=c++11 -qsmp=omp -fopenmp

clean:
	rm -rf ./prog
