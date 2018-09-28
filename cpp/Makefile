CXX = g++
CXXFLAGS = -std=c++11 -fopenmp
LIB = -lbz2 -lpthread
INC=

######
# Main operations
######

all:    iid non_iid restart conditioning

clean:
	rm -f ea_iid ea_non_iid ea_restart ea_conditioning

iid: iid_main.o
iid_main.o: iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) iid_main.cpp $(LIB) -o ea_iid
	
non_iid: non_iid_main.o	
non_iid_main.o: non_iid_main.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) non_iid_main.cpp -o ea_non_iid

restart: restart_main.o
restart_main.o: restart_main.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) restart_main.cpp -o ea_restart

conditioning: conditioning_main.o
conditioning_main.o: conditioning_main.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) conditioning_main.cpp -o ea_conditioning

run_iid:
	./ea_iid ../bin/truerand_8bit.bin 8 -i -t -v

run_non_iid:
	./ea_non_iid ../bin/truerand_8bit.bin 8 -i -t -v
