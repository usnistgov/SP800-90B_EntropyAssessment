#pragma once

#include "utils.h"

/*
* ---------------------------------------------
*		LONGEST REPEATED SUBSTRING HELPERS
* ---------------------------------------------
*/

void find_substrings(const byte text[], int substr_len, map<vector<byte>, vector<int>> &indexes, const int sample_size){

	// First iteration
	if(substr_len == 2){

		// Store all 2-tuples that appear in the text
		for(int i = 0; i < sample_size-1; i++){
			indexes[substr(text, i, substr_len, sample_size)].push_back(i);
		}

	// All other iterations, don't just find all n-tuples naively like 2-tuples
	// Any (n+1)-tuple must build upon an n-tuple, so just take the existing indexes and build on them
	}else{

		vector<int> good_indexes;
		map<vector<byte>, vector<int>>::iterator itr;

		// Store all the indexes of n-tuple substrings that occur more than once
		for(itr = indexes.begin(); itr != indexes.end(); ++itr){
			for(unsigned int i = 0; i < itr->second.size(); i++){
				good_indexes.push_back(itr->second[i]);
			}

		}

		indexes.clear();

		// Extend the n-tuples to (n+1)-tuples and store like before
		for(unsigned int i = 0; i < good_indexes.size(); i++){
			indexes[substr(text, good_indexes[i], substr_len, sample_size)].push_back(good_indexes[i]);
		}
	}
}

void erase_substrings(map<vector<byte>, vector<int>> &indexes){

	// Prune the map of any substrings that occur only once
	map<vector<byte>, vector<int>>::iterator itr = indexes.begin();
	while(itr != indexes.end()){
		if(itr->second.size() < 2){
			indexes.erase(itr++);
		}else{
			++itr;
		}
	}
}

int len_LRS(const byte text[], const int remaining, const int substr_len, const int sample_size){

	// String is the substring we are looking at, vector stores the indexes those substrings begin at
	map<vector<byte>, vector<int>> indexes;
	int len = substr_len;
	int index_size = remaining + 1;

	// Progressively grow the length of the n-tuples to look for
	while(index_size >= remaining){
		len++;
		find_substrings(text, len, indexes, sample_size);
		erase_substrings(indexes);
		index_size = indexes.size();
	}

	if(remaining != 1){ return len-2; }

	// We advance a bit further than we need to
	return len;
}

void count_tuples(const byte data[], const int length, map<vector<byte>, int> &tuples, const int sample_size){

	for(int i = 0; i < sample_size-length; i++){
		vector<byte> substring = substr(data, i, length, sample_size);
		if(tuples.find(substring) == tuples.end()){
			tuples[substring] = 1;
		}else{
			tuples[substring]++;
		}
	}
}

/*
* ---------------------------------------------
*			 HELPER FUNCTIONS
* ---------------------------------------------
*/

void calc_collision_proportion(const vector<double> &p, double &p_col){
	
	for(unsigned int i = 0; i < p.size(); i++){
		p_col += pow(p[i], 2);
	}
}

/*
* ---------------------------------------------
* 		  			 TEST
* ---------------------------------------------
*/

bool len_LRS_test(const byte data[], const int sample_size, const int alphabet_size, const bool verbose){

	vector<double> p(alphabet_size, 0.0);
	calc_proportions(data, p, sample_size);

	double p_col = 0.0;
	calc_collision_proportion(p, p_col);

	// Calculate the number of overlapping substrings of the same length as the longest repeated substring
	int lrs = len_LRS(data, 1, 1, sample_size);
	int n = sample_size - lrs + 1;
	long int overlap = n_choose_2(n);

	double pr_x = 1 - pow(1 - pow(p_col, lrs), overlap);

	if(verbose){
		cout << "Longest Repeated Substring results" << endl;
		cout << "\tP_col: " << p_col << endl;
		cout << "\tLength of LRS: " << lrs << endl;
		cout << "\tPr(X >= 1): " << pr_x << endl;
	}

	return (pr_x >= 0.001);
}

double LRS_test_noniid(const byte data[], const int sample_size, const int alphabet_size, const long u){

	// int u = len_LRS(data, 20, 1, sample_size);
	int v = len_LRS(data, 1, u, sample_size);

	if(v < u){
		cout << "Error in LRS. Aborting." << endl;
		return 0.0;
	}

	vector<double> p;
	for(int i = (int)u; i <= v; i++){
		map<vector<byte>, int> tuples;
		count_tuples(data, i, tuples, sample_size);

		map<vector<byte>, int>::iterator itr;
		int numer = 0;
		for(itr = tuples.begin(); itr != tuples.end(); ++itr){
			if(itr->second != 1){
				numer += n_choose_2(itr->second);
			}
		}

		long int denom = n_choose_2(sample_size-i+1);
		p.push_back(pow(numer/(double)denom, 1.0/i));
	}

	double p_max = 0.0;
	for(unsigned int i = 0; i < p.size(); i++){
		if(p[i] > p_max){
			p_max = p[i];
		}
	}

	return -log2(p_max);
}
