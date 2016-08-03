#include "../shared/utils.h"

double multi_mmc_test(const byte data[]){

	// Step 1
	int D = 16;
	const int N = SIZE - 2;
	array<int, 16> subpredict;
	vector<int> correct(N, false);

	// Step 2
	//  d        x                    y     M_d[x, y]
	map<int, map<array<byte, 16>, map<byte, int>>> M;

	// Step 3
	array<int, 16> scoreboard = {0};
	int winner = 0;

	// Step 4
	for(int i = 3; i < SIZE+1; i++){

		array<byte, 16> substring;

		#ifdef VERBOSE
		if(i % 10000 == 0){
			cout << "\rMultiMMC Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
		}
		#endif

		// Step 4a
		for(int d = 0; d < D; d++){
			substring = fast_substr(data, i-d-1, d);

			if(d+1 < i-1){
				if(M[d].find(substring) == M[d].end()){
					M[d][substring][data[i-1]] = 0;
				}

				M[d][substring][data[i-1]]++;
			}
		}

		// Step 4b
		for(int d = 0; d < D; d++){
			int ymax = 0;
			substring = fast_substr(data, i-d, d);

			if(M[d].find(substring) != M[d].end()){

				byte key = max_map(M[d][substring]);
				ymax = M[d][substring][key];

				if(M[d][substring][ymax] > 0){
					subpredict[d] = ymax;
				}
			}else{
				subpredict[d] = -1;
			}
		}

		// Step 4c
		byte predict = subpredict[winner];

		// Step 4d
		if(predict == data[i-1]){
			correct[i-3] = true;
		}

		// Step 4e
		for(int d = 0; d < D; d++){
			if(subpredict[d] == data[i-1]){
				scoreboard[d]++;

				if(scoreboard[d] >= scoreboard[winner]){
					winner = d;
				}
			}
		}
	}

	// Offsets from flush printing to place carrige on next line
	#ifdef VERBOSE
	cout << endl;
	#endif

	// Step 5
	int C = sum(correct);

	// Step 6
	double p_avg = calc_p_avg(C, N);

	// Step 7
	double p_run = calc_run(correct);

	#ifdef VERBOSE
	cout << "Correct: " << C << endl;
	cout << "P_avg (global): " << p_avg << endl;
	cout << "P_run (local): " << p_run << endl;
	#endif

	// Step 8
	return -log2(max(p_avg, p_run));
}