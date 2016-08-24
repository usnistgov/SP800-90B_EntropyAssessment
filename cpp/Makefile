iid:
	g++ -std=c++11 iid_main.cpp -lbz2 -pthread

fast-iid:
	g++ -std=c++11 -O2 iid_main.cpp -lbz2 -pthread

run-iid:
	./a.out ../bin/truerand_8bit.bin 8 4 -v

non-iid:
	g++ -std=c++11 non_iid_main.cpp

fast-non-iid:
	g++ -std=c++11 -O2 non_iid_main.cpp

run-non-iid:
	./a.out ../bin/truerand_8bit.bin 8 -v
