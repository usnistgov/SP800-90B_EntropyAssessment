#include "../shared/utils.h"

// Much faster method for the first go through (finds the frequency of the mode)
int count_tuple(const byte data[]){

	array<int, 256> freq = {0};
	int mode_freq = -1;
	for(int i = 0; i < SIZE; i++){
		freq[data[i]]++;
		if(freq[data[i]] > mode_freq){
			mode_freq = freq[data[i]];
		}
	}

	return mode_freq;
}

// Slower method that counts frequency of tuples
int count_tuples(const byte data[], const int tuple_size){

	map<vector<byte>, int> count;

	for(int i = 0; i < SIZE+1-tuple_size; i++){

		vector<byte> tuple;
		for(int j = 0; j < tuple_size; j++){
			tuple.push_back(data[i+j]);
		}

		if(count.find(tuple) == count.end()){
			count[tuple] = 1;
		}else{
			count[tuple]++;
		}
	}

	int max = 0;
	map<vector<byte>, int>::iterator itr;
	for(itr = count.begin(); itr != count.end(); ++itr){
		if(itr->second > max){
			max = itr->second;
		}
	}

	return max;
}

double t_tuple_test(const byte data[]){
	
	vector<int> Q;
	int num_tuples = SIZE;
	int tuple_size = 1;

	while(true){

		if(tuple_size == 1){
			num_tuples = count_tuple(data);
		}else{
			num_tuples = count_tuples(data, tuple_size);
		}

		if(num_tuples >= 35){
			Q.push_back(num_tuples);
			tuple_size++;
		}else{
			break;
		}
	}

	vector<double> P, P_max;
	for(int i = 0; i < Q.size(); i++){
		P.push_back(divide(Q[i], SIZE-i));
		P_max.push_back(pow(P[i], divide(1, i+1)));
	}

	return -log2(max_vector(P_max));
}