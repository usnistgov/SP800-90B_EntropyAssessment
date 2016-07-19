#include <stdio.h>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <time.h>
#include <stdlib.h>
using namespace std;

typedef unsigned char byte;

#define SIZE 1000

class Node{
public:
    int by;
    unordered_map<int, Node*> children;
    int indexes;
    Node(int b):by(b){indexes = 0;}
    ~Node(){
        unordered_map<int, Node*>::iterator itr;
        for(itr = children.begin(); itr != children.end(); ++itr)
            delete(itr->second);
    }
};

int max_len = 0;

void insert_in_suffix_tree(Node* root, list<byte> str, int level=0){
    
    root->indexes++;

    if(root->indexes > 1 && max_len < level)
    {
        max_len = level;
    }

    if(str.empty()) return;

    Node* child;
    if(root->children.count(str.front()) == 0) {
        child = new Node(str.front());
        root->children[str.front()] = child;
    } else {
        child = root->children[str.front()];
    }

    str.pop_front();
    insert_in_suffix_tree(child, str, level+1);
}

int main(){

    Node* root = new Node(0);
    list<byte> byte_str;

    srand(time(NULL));
    for(int i = 0; i < SIZE; i++){
        byte_str.push_back(rand() % 256);
    }

    byte_str.push_back(-1);

    for(int i = 0; i < SIZE; i++){
        byte_str.pop_front();
        insert_in_suffix_tree(root, byte_str);
    }
    
    delete(root);

    printf("LRS: %d\n", max_len);

    return 1;
}