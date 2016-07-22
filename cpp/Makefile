all: compile run

compile:
	g++ -std=c++11 iid_main.cpp -lbz2

verbose: compile-verbose run

compile-verbose:
	g++ -std=c++11 iid_main.cpp -lbz2 -DVERBOSE

run:
	./a.out
