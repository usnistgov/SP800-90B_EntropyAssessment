#pragma once

#include "utils.h"

#include <unordered_map>

/*
* ---------------------------------------------
*			 SUFFIX TREE DEFINITION
* ---------------------------------------------
*/

int max_len = 0;

class Node{
	public:
		byte b;
		unordered_map<byte, Node*> children;
		vector<int> indexes;
		Node(byte b_in):b(b_in){}
		~Node(){
			unordered_map<byte, Node*>::iterator itr;
			for(itr = children.begin(); itr != children.end(); ++itr)
				delete(itr->second);
		}
};

void insert_in_suffix_tree(Node* root, vector<byte> byte_str, int index, int level = 0){

	root->indexes.push_back(index);

	if(root->indexes.size() > 1 && max_len < level){
		max_len = level;
	}

	if(byte_str.empty()) return;

	Node *child;
	if(root->children.count(byte_str[0]) == 0){
		child = new Node(byte_str[0]);
		root->children[byte_str[0]] = child;
	}else{
		child = root->children[byte_str[0]];
	}

	vector<byte>::iterator front_itr = byte_str.begin()+1, back_itr = byte_str.end();
	byte_str = vector<byte>(front_itr, back_itr);
	insert_in_suffix_tree(child, byte_str, index, level+1);
}

int len_LRS(const byte data[]){

	Node* root = new Node(0);
	vector<byte> byte_str;
	
	int factor = 250;

	for(long int i = 0; i < SIZE/factor; i++){
		byte_str.push_back(data[i]);
	}

	cout << "Vector built" << endl;

	vector<byte>::iterator front_itr, back_itr = byte_str.end();
	vector<byte> next_byte_str;
	for(int i = 0; i < SIZE/factor; i++){

	 	front_itr = byte_str.begin() + i;

	 	next_byte_str = vector<byte>(front_itr, back_itr);
	 	insert_in_suffix_tree(root, next_byte_str, i);
	}
	
	delete(root);

	return max_len;
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

bool len_LRS_test(const byte data[]){

	vector<double> p(256, 0.0);
	calc_proportions(data, p); // borrowed from chi_square_tests.h

	double p_col = 0.0;
	calc_collision_proportion(p, p_col);

	cout << len_LRS(data) << endl;

	return true;
}