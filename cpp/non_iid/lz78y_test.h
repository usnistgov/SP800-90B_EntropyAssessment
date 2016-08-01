#include "../shared/utils.h"

double LZ78Y_test(const byte data[]){

	// Step 1
	int B = 16;
	int N = SIZE - B - 1;
	vector<int> correct(SIZE-B-1, false);	
	int max_dict_size = 65536;

	// Step 2
	map<array<byte, 16>, map<byte, int>> D;
	int D_size = 0;

	// Step 3
	for(int i = B+1; i < SIZE; i++){

		// Pull entire region for faster look ups
		array<byte, 16> substring, subregion = fast_substr(data, i-B-1, B);

		#ifdef VERBOSE
			if(i % 10000 == 0){
				cout << "\rLZ78Y Test: " << (i/(double)SIZE)*100 << "% complete" << flush;
			}
		#endif

		// Step 3a
		for(int j = B; j >= 1; j--){

			substring = fast_substr(subregion.data(), B-j, j);

			if(D_size < max_dict_size){
				if(D.find(substring) == D.end()){
					D[substring][data[i-1]]++;
				}else{
					D[substring][data[i-1]] = 0;
					D_size++;
				}
			}
		}

		// Step 3b
		int max_count = 0;
		byte prediction;

		for(int j = B; j >= 1; j--){

			substring = fast_substr(subregion.data(), B-j+1, j);

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

	#ifdef VERBOSE
		cout << endl;
	#endif

	// Step 4
	int C = sum(correct);
	double p_avg = calc_p_avg(C, N);
	
	// Step 5
	double p_run = calc_run(correct);

	cout << "Correct: " << C << endl;
	cout << "P_avg (global): " << p_avg << endl;
	cout << "P_run (local): " << p_run << endl;

	// Step 6
	return -log2(max(p_avg, p_run));
}