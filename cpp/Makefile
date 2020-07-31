#CXX = g++
#CXXFLAGS = -std=c++11 -fopenmp -O2 -msse2 -ffloat-store -march=native
CXX = clang++
# -static-libstdc++ -static-libgcc
CXXFLAGS = -Wno-padded -Wno-disabled-macro-expansion -Wno-gnu-statement-expression -Wno-bad-function-cast -fopenmp -O1 -fsanitize=address -fsanitize=undefined -fdenormal-fp-math=ieee -msse2 -march=native -fforce-enable-int128 -L C:\Users\skb1\csdprojects\libdivsufsort-master\out\build\x86-Debug\lib -L C:\Users\skb1\Downloads\getopt_mb_uni_vc10\getopt_mb_uni_vc10\bin\debug\dll -L C:\Users\skb1\csdprojects\bzip2-1.0.6-master_secondcopy -pthread -Xlinker /NODEFAULTLIB:msvcrt

# -Xlinker /NODEFAULTLIB:libcmt -Xlinker /NODEFAULTLIB:ucrt

# -Xlinker /NODEFAULTLIB
# 
# bzip was compiled using nmake and makefile.msc using the Visual Studio command prompt
#static analysis in clang using
#scan-build-8 --use-c++=/usr/bin/clang++-8 make
LIB = -lgetopt -llibbz2 -ldivsufsort 
# -lmsvcrt -lucrt -lmsvcprt -lvcruntime
INC=-I C:\Users\skb1\Downloads\getopt_mb_uni_vc10\getopt_mb_uni_vc10\getopt_mb_uni_vc10_dll -I C:\Users\skb1\csdprojects\unistd_h-master\include -I C:\Users\skb1\csdprojects\libdivsufsort-master\out\build\x86-Debug\include -I C:\Users\skb1\csdprojects\bzip2-1.0.6-master_secondcopy


######
# Main operations
######

all:    iid non_iid restart conditioning transpose

clean:
	rm -f ea_iid ea_non_iid ea_restart ea_conditioning ea_transpose selftest/*.res

iid: iid_main.o
iid_main.o: iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) iid_main.cpp -o ea_iid.exe $(LIB)
	
non_iid: non_iid_main.o	
non_iid_main.o: non_iid_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) non_iid_main.cpp -o ea_non_iid.exe $(LIB)

restart: restart_main.o
restart_main.o: restart_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) restart_main.cpp -o ea_restart.exe $(LIB)

conditioning: conditioning_main.o
conditioning_main.o: conditioning_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) conditioning_main.cpp -o ea_conditioning.exe $(LIB)

transpose: transpose_main.o
transpose_main.o: transpose_main.cpp
	$(CXX) $(CXXFLAGS) $(INC) transpose_main.cpp -o ea_transpose.exe $(LIB)
