#include "../shared/utils.h"

double LZ78Y_test(const byte data[], const bool verbose){

	// Step 1
	int B = 16;
	int N = SIZE - B - 1;
	vector<int> correct(SIZE-B-1, false);	
	int max_dict_size = 65536;

	// Step 2
	map<array<byte, 16>, map<byte, int>> D;
	int D_size = 0;

	// Step 3
	array<byte, 16> substring;
	for(int i = B+1; i < SIZE; i++){

		if(verbose){
			if(i % 10000 == 0){
				cout << "\rLZ78Y Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
			}
		}

		// Step 3a
		for(int j = B; j >= 1; j--){
			substring = fast_substr(data, i-j-1, j);

			if(D_size < max_dict_size){
				if(D.find(substring) == D.end()){
					D[substring][data[i-1]]++;
				}else{
					D[substring][data[i-1]] = 1;
					D_size++;
				}
			}
		}

		// Step 3b
		int max_count = 0;
		byte prediction = 0;

		for(int j = B; j >= 1; j--){
			substring = fast_substr(data, i-j, j);

			if(D.find(substring) != D.end()){
				byte y = max_map(D[substring]);
				if(D[substring][y] > max_count){
					prediction = y;
					max_count = D[substring][y];
				}
			}
		}

		// Step 3c
		if(prediction == data[i]){
			correct[i-B-1] = true;
		}
	}

	if(verbose){
		cout << endl;
	}

	// Step 4
	int C = sum(correct);
	double p_avg = calc_p_avg(C, N);
	
	// Step 5
	double p_run = calc_run(correct);

	if(verbose){
		cout << "\tCorrect: " << C << endl;
		cout << "\tP_avg (global): " << p_avg << endl;
		cout << "\tP_run (local): " << p_run << endl;
	}

	// Step 6
	return -log2(max(p_avg, p_run));
}