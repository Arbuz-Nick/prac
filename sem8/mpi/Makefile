all: main

problem: *.cpp *.h
	mpicxx problem_main.cpp matrix.cpp -o prog -std=c++11

main: *.cpp *.h
	mpicxx main.cpp matrix.cpp -o prog -std=c++11

clean:
	rm -rf ./prog
