CXX = g++
CXXFLAGS = -std=c++11 -fopenmp
CDDFLAGS = -std=c++11 -ggdb -fopenmp
TESTFLAGS = -ftree-vectorize -msse2 -Ofast

LIB = -lbz2 -lpthread
INC=-I/usr/include -I/usr/include/c++/4.9 -I/usr/include/c++/4.9/backward \
 -I/usr/include/x86_64-linux-gnu -I/usr/include/x86_64-linux-gnu/c++/4.9 \
 -I/usr/lib/gcc/x86_64-linux-gnu/4.9/include -I/usr/lib/gcc/x86_64-linux-gnu/4.9/include-fixed \
 -I/usr/local/include

######
# Main operations
######

iid: iid_main.o
iid_main.o: iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) iid_main.cpp $(LIB) -o ea_iid
	
non_iid: non_iid_main.o	
non_iid_main.o: non_iid_main.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) non_iid_main.cpp -o ea_non_iid

restart: restart_main.o
restart_main.o: restart_main.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) restart_main.cpp -o ea_restart

run_iid:
	./ea_iid ../bin/truerand_8bit.bin 8 -i -t -v

run_non_iid:
	./ea_non_iid ../bin/truerand_8bit.bin 8 -i -t -v


######
# Test operations
######

chi_test: chi_test_main.o
chi_test_main.o: test/chi_test_main.cpp
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) $(INC) test/chi_test_main.cpp $(LIB) -o chi_test

run_chi_test:
	./chi_test

lrs_test: lrs_test_main.o
lrs_test_main.o: test/lrs_test_main.cpp
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) $(INC) test/lrs_test_main.cpp $(LIB) -o lrs_test

run_lrs_test:
	./lrs_test

most_common_test: most_common_test_main.o
most_common_test_main.o: test/most_common_test_main.cpp
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) $(INC) test/most_common_test_main.cpp $(LIB) -o most_common_test

run_most_common_test:
	./most_common_test

permutation_test: permutation_test_main.o
permutation_test_main.o: test/permutation_test_main.cpp
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) $(INC) test/permutation_test_main.cpp $(LIB) -o permutation_test

run_permutation_test:
	./permutation_test ../bin/truerand_8bit.bin 8 1000000 1 -v
