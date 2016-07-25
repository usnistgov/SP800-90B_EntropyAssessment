iid: compile-iid run

compile-iid:
	g++ -std=c++11 iid_main.cpp -lbz2

verbose-iid: compile-verbose-iid run

compile-verbose-iid:
	g++ -std=c++11 iid_main.cpp -lbz2 -DVERBOSE

run:
	./a.out

non-iid: compile-non-iid run

compile-non-iid:
	g++ -std=c++11 non_iid_main.cpp

verbose-non-iid: compile-verbose-non-iid run

compile-verbose-non-iid:
	g++ -std=c++11 non_iid_main.cpp -DVERBOSE
