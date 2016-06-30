#include <time.h>
#include <stdlib.h>
#include <iostream>

#define SIZE 1000000
#define PERMS 100

#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)

int main(){

	int vals[SIZE];

	for(int i = 0; i < SIZE; ++i){
		vals[i] = i%256;
	}

	srand(time(NULL));
	for(int j = 0; j < PERMS; ++j){
		int r;
		for(int i = SIZE-1; i > 0; --i){
			r = (rand() / (float)RAND_MAX) * (i+1);
			SWAP(vals[r], vals[i]);
		}
	}

	for(int i = 0; i < 10; i++){
		std::cout << rand() << std::endl;
	}
}