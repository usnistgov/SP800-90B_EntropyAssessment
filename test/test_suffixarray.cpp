#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cstring>
#include <iomanip>

using namespace std;

//typedef unsigned char byte;

#define SIZE 1000000

class SuffixArray{
	private:
		string *text;
		int length;
		int *index;
		string *suffix;

	public:
		SuffixArray(string text){
			
			this->text = new string[text.length()];

			for(int i = 0; i < text.length(); i++){
				this->text[i] = text.substr(i, 1);
			}

			this->length = text.length();
			this->index = new int[length];

			for(int i = 0; i < length; i++){
				index[i] = i;
			}

			suffix = new string[length];
		}

		void createSuffixArray(){

			for(int i = 0; i < length; i++){

				string text = "";

				for(int j = i; j < length; j++){
					text += this->text[j];
				}

				suffix[i] = text;
			}

			int back;
			for(int i = 1; i < length; i++){
				string key = suffix[i];
				int key_index = index[i];

				for(back = i-1; back >= 0; back--){

					if(suffix[back].compare(key) > 0){

						suffix[back+1] = suffix[back];
						index[back+1] = index[back];
					}else{
						break;
					}
				}

				suffix[back+1] = key;
				index[back+1] = key_index;
			}
		}

		void print(){
			int longest = 0;

			for(int i = 0; i < length; i++){
				if(suffix[i].length() > longest){
					longest = suffix[i].length();
				}
			}

			cout << setw(longest+2) << "SUFFIX";
			cout << " INDEX" << endl;
			for(int i = 0; i < length; i++){
				cout << setw(longest+2) << suffix[i] << " ";
				cout << index[i] << endl;
			}
		}

		void longestCommonPrefix(){

			int ret = 0;
			string val = "";

			// take each suffix and see if it exists letter by letter in the next suffix
			string cur, next;
			for(int i = 0; i < length-1; i++){
				cur = suffix[i];
				next = suffix[i+1];

				for(int j = 0; j < min(cur.length(), next.length()); j++){
					if(cur[j] == next[j]){
						if(ret < (j+1)){
							ret = j+1;
							val = cur.substr(0, j+1);
						}
					}else{
						break;
					}
				}
			}

			cout << ret << " from \"" << val << "\"" << endl;
		}
};

int main(){

	string text = "mississippi";
	SuffixArray sa(text);
	sa.createSuffixArray();
	sa.print();
	sa.longestCommonPrefix();

	return 0;
}