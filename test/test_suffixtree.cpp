#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
using namespace std;

typedef unsigned char byte;

class Node
{
public:
    byte ch;
    unordered_map<byte, Node*> children;
    vector<int> indexes; //store the indexes of the substring from where it starts
    Node(byte c):ch(c){}
};

int maxLen = 0;
string maxStr = "";

void insertInSuffixTree(Node* root, vector<byte> str, int index, vector<byte> originalSuffix, int level=0)
{
    root->indexes.push_back(index);

    // it is repeated and length is greater than maxLen
    // then store the substring
    if(root->indexes.size() > 1 && maxLen < level)
    {
        maxLen = level;
        //maxStr = originalSuffix.substr(0, level);
    }

    if(str.empty()) return;

    Node* child;
    if(root->children.count(str[0]) == 0) {
        child = new Node(str[0]);
        root->children[str[0]] = child;
    } else {
        child = root->children[str[0]];
    }

    vector<byte>::iterator itr = str.begin();
    str.erase(itr);
    insertInSuffixTree(child, str, index, originalSuffix, level+1);
}

int main()
{
    byte str[] = {0, 1, 2, 3, 4, 2, 3, 4, 2, 1};
    Node* root = new Node(0);

    vector<byte> vb;
    for(int i = 0; i < 10; i++){
        vb.push_back(str[i]);
    }

    //insert all substring in suffix tree
    for(int i=0; i<10; i++){
        vector<byte> vb2 = vb;
        for(int j = 0; j < i; j++){
            vector<byte>::iterator itr = vb2.begin();
            vb2.erase(itr);
        }
        insertInSuffixTree(root, vb2, i, vb2);
    }

    cout << maxLen << endl;

    return 1;
}