#include "../shared/utils.h"

double lag_test(const byte data[]){
	
	// Step 1
	int D = 128;
	int N = SIZE-1;
	vector<int> correct(N, false);
	array<int, 128> lag;

	// Step 2
	array<int, 128> scoreboard;
	int winner = 0;

	// Step 3
	for(int i = 2; i < SIZE+1; i++){

		#ifdef VERBOSE
		if(i % 10000 == 0){
			cout << "\rLag Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
		}
		#endif

		// Step 3a
		for(int d = 1; d < D+1; d++){
			if(d < i){
				lag[d-1] = data[i-d-1];
			}else{
				lag[d-1] = -1;
			}
		}

		// Step 3b-c
		correct[i-2] = (lag[winner] == data[i-1]);

		// Step 3d
		for(int d = 0; d < D; d++){
			if(lag[d] == data[i-1]){
				scoreboard[d]++;
				if(scoreboard[d] >= scoreboard[winner]){
					winner = d;
				}
			}
		}
	}

	#ifdef VERBOSE
	cout << endl;
	#endif

	// Step 4
	int C = sum(correct);

	// Step 5
	double p_avg = calc_p_avg(C, N);

	// Step 6
	double p_run = calc_run(correct);

	#ifdef VERBOSE
	cout << "\tCorrect: " << C << endl;
	cout << "\tP_avg (global): " << p_avg << endl;
	cout << "\tP_run (local): " << p_run << endl;
	#endif

	// Step 7
	return -log2(max(p_avg, p_run));
}