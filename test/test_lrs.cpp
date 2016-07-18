#include <iostream>
#include <string>

using namespace std;

int find_LRS(string str){

	int n = str.length();

	int dp[n+1][n+1];

	for(int i = 0; i < n+1; i++){
		for(int j = 0; j < n+1; j++){
			dp[i][j] = 0;
		}
	}

	for(int i = 1; i < n+1; i++){
		for(int j = 1; j < n+1; j++){
			if(str[i-1] == str[j-1] && i != j){
				dp[i][j] = 1 + dp[i-1][j-1];
			}else{
				dp[i][j] = max(dp[i][j-1], dp[i-1][j]);
			}
		}
	}

	for(int i = 0; i < n+1; i++){
		for(int j = 0; j < n+1; j++){
			cout << dp[i][j] << " ";
		}
		cout << endl;
	}

	return dp[n][n];
}

int main(){
	string str = "11231234";
	cout << find_LRS(str) << endl;
	return 0;
}