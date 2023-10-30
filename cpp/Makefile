ARCH ?= x86

CXX ?= $(CROSS_COMPILE)g++

CXXFLAGS = -std=c++11 -fopenmp -O2 -ffloat-store -I/usr/include/jsoncpp
ifeq ($(ARCH),x86)
CXXFLAGS += -msse2 -march=native
endif

#CXX = clang++-15
#CXXFLAGS = -g -Wno-padded -Wno-disabled-macro-expansion -Wno-gnu-statement-expression -Wno-bad-function-cast -fopenmp -O1 -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer -fdenormal-fp-math=ieee -msse2 -march=native -I/usr/include/jsoncpp
#static analysis in clang using
#scan-build-15 --use-c++=/usr/bin/clang++-15 make
LIB = -lbz2 -lpthread -ldivsufsort -ldivsufsort64
COND_LIB = -lmpfr -lgmp
SHARED_LIB = -ljsoncpp -lcrypto
INC =

######
# Main operations
######

all:    iid non_iid restart conditioning transpose

clean:
	rm -f ea_iid ea_non_iid ea_restart ea_conditioning ea_transpose selftest/*.res

iid: iid_main.o
iid_main.o: iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) iid_main.cpp -o ea_iid $(LIB) $(SHARED_LIB)
	
non_iid: non_iid_main.o	
non_iid_main.o: non_iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) non_iid_main.cpp -o ea_non_iid $(LIB) $(SHARED_LIB)

restart: restart_main.o
restart_main.o: restart_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) restart_main.cpp -o ea_restart $(LIB) $(SHARED_LIB)

conditioning: conditioning_main.o
conditioning_main.o: conditioning_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) conditioning_main.cpp -o ea_conditioning $(LIB) $(COND_LIB) $(SHARED_LIB)

transpose: transpose_main.o
transpose_main.o: transpose_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) transpose_main.cpp -o ea_transpose $(LIB) $(SHARED_LIB)
