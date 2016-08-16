#include "../shared/utils.h"

byte most_common_recent(const byte data[], const int pos, const int length){

	byte current_mode = 0;
	int mode_frequency = -1;
	array<int, 256> freq = {0};

	for(int i = 0; i < length; i++){
		byte cur_val = data[pos+i];
		freq[cur_val]++;

		if(freq[cur_val] >= mode_frequency){
			current_mode = cur_val;
			mode_frequency = freq[cur_val];
		}
	}

	return current_mode;
}

double multi_mcw_test(const byte data[], const bool verbose){
	
	// Step 1
	array<int, 4> w = {63, 255, 1023, 4095};
	int N = SIZE - w[0];
	vector<int> correct(N, false);

	// Step 2
	array<int, 4> scoreboard = {0};
	array<int, 4> frequent = {-1, -1, -1, -1};
	int winner = 0;

	// Step 3
	for(int i = w[0]+1; i < SIZE+1; i++){

		if(verbose){
			if(i % 10000 == 0){
				cout << "\rMulti Most Common in Window (MultiMCW) Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
			}
		}

		// Step 3a
		for(int j = 0; j < 4; j++){
			if(i > w[j]+1){
				frequent[j] = most_common_recent(data, i-w[j]-1, w[j]);
			}else{
				frequent[j] = -1;
			}
		}

		// Step 3b-c
		correct[i-w[0]-1] = (frequent[winner] == data[i-1]);

		// Step 3d
		for(int j = 0; j < 4; j++){
			if(frequent[j] == data[i-1]){
				scoreboard[j]++;
				if(scoreboard[j] >= scoreboard[winner]){
					winner = j;
				}
			}
		}
	}

	if(verbose){
		cout << endl;
	}

	// Step 4
	int C = sum(correct);

	// Step 5
	double p_avg = calc_p_avg(C, N);

	// Step 6
	double p_run = calc_run(correct);

	if(verbose){
		cout << "\tCorrect: " << C << endl;
		cout << "\tP_avg (global): " << p_avg << endl;
		cout << "\tP_run (local): " << p_run << endl;
	}

	// Step 7
	return -log2(max(p_avg, p_run));
}