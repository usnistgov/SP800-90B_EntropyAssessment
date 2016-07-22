#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

typedef unsigned char byte;

using namespace std;

vector<byte> substr(byte text[], int pos, int len){

	int size = 10;

	if(pos+len > size){
		len = size - pos;
	}

	vector<byte> substring;

	for(int i = 0; i < len; i++){
		substring.push_back(text[pos+i]);
	}

	return substring;
}

void find_substrings(byte text[], int substr_len, map<vector<byte>, vector<int>> &indexes){
	
	int size = 10;

	// First iteration
	if(substr_len == 2){

		// Store all 2-tuples that appear in the text
		for(int i = 0; i < size-1; i++){
			indexes[substr(text, i, substr_len)].push_back(i);
		}

	// All other iterations, don't just find all n-tuples naively like 2-tuples
	// Any (n+1)-tuple must build upon an n-tuple, so just take the existing indexes and build on them
	}else{

		vector<int> good_indexes;
		map<vector<byte>, vector<int>>::iterator itr;

		// Store all the indexes of n-tuple substrings that occur more than once
		for(itr = indexes.begin(); itr != indexes.end(); ++itr){
			for(int i = 0; i < itr->second.size(); i++){
				good_indexes.push_back(itr->second[i]);
			}

		}

		indexes.clear();

		// Extend the n-tuples to (n+1)-tuples and store like before
		for(int i = 0; i < good_indexes.size(); i++){
			indexes[substr(text, good_indexes[i], substr_len)].push_back(good_indexes[i]);
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

int len_LRS(byte text[]){

	// String is the substring we are looking at, vector stores the indexes those substrings begin at
	map<vector<byte>, vector<int>> indexes;
	int substr_len = 2;
	
	// Progressively grow the length of the n-tuples to look for
	while(true){
		find_substrings(text, substr_len, indexes);
		erase_substrings(indexes);

		if(indexes.empty()) break;
		substr_len++;
	}

	// We always advance one further than we need to
	return substr_len-1;
}

int main(){

	//ifstream file;
	//string text;
	//string line;

	//file.open("moby_dick.txt");
	//while(getline(file, line)) text += line;
	//file.close();

	//string text = "abcdeabcdeabc";

	byte text[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 3};

	cout << len_LRS(text) << endl;

	return 0;
}