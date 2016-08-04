#include "../shared/utils.h"

byte most_common_recent(const vector<byte> &data){

	map<byte, int> proportions;
	map_init(proportions);

	for(int i = 0; i < data.size(); i++){
		proportions[data[i]]++;
	}

	// Detemine mode(s) of dataset
	set<byte> modes;
	int max = 0;
	map<byte, int>::iterator itr;
	for(itr = proportions.begin(); itr != proportions.end(); ++itr){
		if(itr->second >= max){
			if(itr->second > max){
				modes.clear();
				max = itr->second;
			}

			modes.insert(itr->first);
		}
	}

	// Return the most recent occurance of the mode
	if(modes.size() >= 2){
		for(int i = data.size()-1; i >= 0; i--){
			if(modes.find(data[i]) != modes.end()){
				return data[i];
			}
		}

	}else{
		return *(modes.begin());
	}
}

double multi_mcw_test(const byte data[]){
	
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

		#ifdef VERBOSE
		if(i % 10000 == 0){
			cout << "\rMulti Most Common in Window (MultiMCW) Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
		}
		#endif

		// Step 3a
		for(int j = 0; j < 4; j++){
			if(i > w[j]+1){
				vector<byte> substring = substr(data, i-w[j]-1, w[j]);
				frequent[j] = most_common_recent(substring);
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
	cout << "Correct: " << C << endl;
	cout << "P_avg (global): " << p_avg << endl;
	cout << "P_run (local): " << p_run << endl;
	#endif

	// Step 7
	return -log2(max(p_avg, p_run));
}