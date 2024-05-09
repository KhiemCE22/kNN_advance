
// Online C++ Compiler - Build, Compile and Run your C++ programs online in your favorite browser

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <list>
#include <queue>

using namespace std;
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }
    void print();
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;

public:
    kDTree(int k = 2);
    ~kDTree();

    const kDTree &operator=(const kDTree &other);
    kDTree(const kDTree &other);


    /*recursion method*/
    void inorderRec(kDTreeNode* node) const;
    void preorderRec(kDTreeNode* node) const;
    void postorderRec(kDTreeNode* node) const;
    int heightRec(kDTreeNode* node) const;
    int leafCountRec(kDTreeNode* node) const;
    int nodeCountRec(kDTreeNode* node) const;

    void insertRec(kDTreeNode*& node, const vector<int> & point, int depth);
    bool searchRec(kDTreeNode* node, const vector<int> &point, int depth);
    void buildTreeRec(vector<vector<int>>& pointList, int start, int end, kDTreeNode*& node, int depth);
    
    // ********************
    kDTreeNode* removeRec(const vector<int> &point, kDTreeNode* & node, int depth);
    void nearestNeighbourRec(const vector<int> &target, kDTreeNode *&best, kDTreeNode* node, int depth);
    void kNearestNeighbour(const vector<int> &target, int k, kDTreeNode* node,vector<kDTreeNode *> &bestList, int depth);
    /*official method*/

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const;   
    int nodeCount() const;   
    int leafCount() const;

    
    void insert(const vector<int> &point);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
    
    // bonus for remove method
    kDTreeNode* getMinNodeRec(kDTreeNode*  node, int split_axis, int depth);
    kDTreeNode* getMinNode(int split_axis);
    ///// vissulize
    void visualize() const {
        if (!root)
            return;

        // **Pre-compute node levels for efficiency**
        vector<vector<kDTreeNode*>> levelNodes(heightRec(root));
        getLevelNodes(root, 0, levelNodes);

        // **Print each level with proper spacing**
        for (int i = 0; i < levelNodes.size(); i++) {
            printRow(levelNodes[i], levelNodes.size(), i);
        }
    }

    void getLevelNodes(kDTreeNode* node, int depth, vector<vector<kDTreeNode*>>& levelNodes) const {
        if (!node)
            return;

        levelNodes[depth].push_back(node);

        getLevelNodes(node->left, depth + 1, levelNodes);
        getLevelNodes(node->right, depth + 1, levelNodes);
    }

    void printRow(const vector<kDTreeNode*>& nodes, int totalLevels, int depth) const {
    int spacing = (totalLevels - depth) * 4; // Adjust spacing for better readability

    cout << setw(spacing);
    bool toggle = true; // Start with left
    for (kDTreeNode* node : nodes) {
        if (node) { // Print node data if not nullptr
            if (toggle)
                cout << "/" << setw(spacing - 3);
            else
                cout << "\\" << setw(spacing - 3);

            // Print node data (assuming data[0])
            cout << "(" << node->data[0] << ", " << node->data[1] << ", " << node->data[2] << ")";
        } else { // Print placeholder for empty node
            cout << setw(spacing - 3); // Adjust spacing for placeholder
        }

        toggle = !toggle;
    }
    cout << endl;

    cout << setw(spacing);
    for (kDTreeNode* node : nodes) {
        if (node)
            cout << setw(spacing - 3); // Adjust spacing for node data
        else
            cout << setw(spacing - 3) << "-"; // Print placeholder symbol
    }
    cout << endl;
    }

};

/*Class kDTreeNode*/
void kDTreeNode :: print() {
        cout << "(";
        for (int i = 0; i < data.size(); ++i) {
            if (i<=2){
                cout << data[i];
                if (i < data.size() - 1) {
                    cout << ",";
                }
            }
            else if (i == data.size() - 1 ){
                if (data.size() > 4)
                    cout<<"...,";
                cout<<data[i];
            }
        }
        cout << ")";
}


/*Classs kDTree*/
kDTree :: kDTree(int k){
    this -> k = k;
    this -> root = NULL;
}

/*bonus function for destructor*/
void deleteTree(kDTreeNode*& root){
    if (!root) return;
    deleteTree(root -> left);
    deleteTree(root -> right);
    delete root;
    root = NULL;
}
kDTree :: ~kDTree(){
    deleteTree(this -> root);
}
void kDTree :: inorderRec(kDTreeNode* node) const {
    if (!node) return ;
    inorderRec(node -> left);
    if (node -> left)
        cout<<" ";
    node -> print();
    if (node -> right)
        cout<<" ";
    inorderRec(node -> right);
}
void kDTree :: inorderTraversal()const{
    return inorderRec(this -> root);
}
void kDTree :: preorderRec(kDTreeNode* node)const{
    if (!node) return;
    node -> print();
    if (node -> left)
        cout<<" ";
    preorderRec(node -> left);
    if (node -> right)
        cout<<" ";
    preorderRec(node -> right);
}
void kDTree :: preorderTraversal()const{
    return preorderRec(this -> root);
}
void kDTree :: postorderRec(kDTreeNode* node)const{
    if (!node) return;
    postorderRec(node -> left);
    if (node -> left)
        cout<<" ";
    postorderRec(node -> right);
    if (node -> right)
        cout<<" ";
    node -> print();
}
void kDTree :: postorderTraversal()const{
    return postorderRec(this -> root);
}

int kDTree :: heightRec(kDTreeNode* node) const{
    if (!node) return 0;
    int hLeft = heightRec(node -> left);
    int hRight = heightRec(node -> right);
    if (!node -> left) return hRight + 1;
    else if (!node -> right) return hLeft + 1;
    else 
        return hLeft > hRight ? hLeft + 1 : hRight + 1;
}
int kDTree ::  height() const{
    return heightRec(this -> root);
}

int kDTree :: nodeCountRec(kDTreeNode* node)const{
    if (!node) return 0;
    return 1 + nodeCountRec(node -> left) + nodeCountRec(node -> right);
}
int kDTree :: nodeCount()const{
    return nodeCountRec(this -> root);
}

int kDTree :: leafCountRec(kDTreeNode* node) const{
    if (!node) return 0;
    if (!node -> left && !node -> right) return 1;
    else return leafCountRec(node -> left) + leafCountRec(node -> right);
}
int kDTree :: leafCount()const{
    return leafCountRec(this -> root);
}

void kDTree::insertRec(kDTreeNode*& node, const vector<int> & point, int depth) {
    // if dim of data is not equal k, do not thing
    if (point.size() != this -> k) return;
    if (!node){
        node = new kDTreeNode(point);
        return;
    }
    int dimension = depth % this ->k;
    if (point[dimension] < node-> data[dimension])
        insertRec(node -> left, point, depth + 1);
    else 
        insertRec(node -> right, point, depth + 1);
}
void kDTree::insert(const vector<int> &point){
    return insertRec(this -> root, point, 0);
}

/*bonus function for search*/

bool equal(const vector<int>& point, vector<int> data ){
    if (data.size() != point.size()) return false;
    for (int i = 0; i < data.size(); i++){
        if (data[i] != point[i]) 
            return false;
    }
    return true;
}

bool kDTree:: searchRec(kDTreeNode* node, const vector<int> &point, int depth){
    
    if (!node) return false;
    if (equal(point, node -> data)) return true;

    int dimension = depth % k;
    if (point[dimension] < node -> data[dimension])
        return searchRec(node -> left, point, depth + 1);
    else 
        return searchRec(node -> right, point, depth + 1);

}

bool kDTree :: search(const vector<int> &point){
    // if dim of data is not equal k, return false
    if (point.size() != this -> k) return false;
    return searchRec(this -> root, point, 0);
}

/*
    use merge sort for build tree
*/


void merge(int left, int middle , int right, vector<vector<int>> &pointList, int dimension){
    // deep copy to store left to middle
    vector<vector<int>> leftList(pointList.begin() + left, pointList.begin() + middle + 1);
    vector<vector<int>> rightList(pointList.begin() + middle + 1, pointList.begin() + right + 1);

    int i = 0, j = 0, k = left;
    int n1 = leftList.size();
    int n2 = rightList.size();

    while (i < n1 && j < n2){
        if (leftList[i][dimension] <= rightList[j][dimension]){
            pointList[k++] = leftList[i++];
            
        }
        else {
            pointList[k++] = rightList[j++];
        }
    }

    while (i < n1){
        pointList[k++] = leftList[i++];
    }
    while (j < n2){
        pointList[k++] = rightList[j++];
    }


}


void mergeSort(vector<vector<int>> & pointList, int start, int end, int dimension){
    if (start < end){
        int middle= start + (end - start)/2;
        mergeSort(pointList,start,middle,dimension);
        mergeSort(pointList,middle+1,end,dimension);
        merge(start, middle, end, pointList, dimension);
    }
    
}

void kDTree :: buildTreeRec( vector<vector<int>> &pointList,int start, int end, kDTreeNode* &node, int depth){
    if (start > end) return;
    int dimension = depth % k;
    mergeSort(pointList, start, end, dimension);
    int median = start + (end - start)/2;
    

    while(median > start && pointList[median][dimension] == pointList[median-1][dimension]){
        median--;
    }
    vector<int> data = pointList[median];
    node = new kDTreeNode(data);

    buildTreeRec(pointList,start,median-1,node -> left, depth + 1);
    buildTreeRec(pointList,median + 1,end,node -> right, depth + 1);    

}
void kDTree :: buildTree(const vector<vector<int>> &pointList){
    // destree 
    vector<vector<int>> points = pointList;
    if (root)
        deleteTree(this -> root);
    buildTreeRec(points, 0, points.size() - 1, this -> root, 0);

}
kDTreeNode* kDTree :: getMinNodeRec(kDTreeNode*  node, int split_axis, int depth){
    if (!node || (!node -> left && !node -> right))
        return node;    
    int dimesion = depth % k;
    if (dimesion == split_axis){
        if (!node -> left)
            return node;
        else 
            return getMinNodeRec(node -> left, split_axis, depth + 1);
    }
    else {
        kDTreeNode* min_left =  getMinNodeRec(node -> left, split_axis, depth + 1);
        kDTreeNode* min_right = getMinNodeRec(node -> right, split_axis, depth + 1);
        if (!min_left)
            return node -> data[split_axis] <= min_right->data[split_axis] ? node : min_right;
        else if (!min_right)
            return node -> data[split_axis] <= min_left->data[split_axis] ? node : min_left;
        else if (node -> data[split_axis] <= min_left->data[split_axis] && node -> data[split_axis] <= min_right->data[split_axis])
            return node;
        else 
            return min_left -> data[split_axis] <= min_right -> data[split_axis] ? min_left : min_right;
    }
}
kDTreeNode* kDTree :: getMinNode(int split_axis){
    return getMinNodeRec(this -> root, split_axis, 0);
}



kDTreeNode* kDTree :: removeRec(const vector<int> &point, kDTreeNode* & node, int depth){
    if (!node) return NULL;
    int dimension =  depth % k;
    if (equal(point, node -> data)){
        // if right subtree not nullptr 
        if (node -> right){
            kDTreeNode* succesor = getMinNodeRec(node -> right, dimension, depth + 1);
            node -> data = succesor -> data;
            node -> right = removeRec(succesor -> data, node -> right, depth + 1);
        }
        // if right subtree is null and left subtree is not null
        else if (node -> left){
            kDTreeNode* sucessor = getMinNodeRec(node -> left, dimension, depth + 1);
            node -> data = sucessor -> data;
            node -> right = removeRec(sucessor -> data, node -> left, depth + 1);
            node -> left = NULL;            
        }
        // if is a leaf
        else {
            delete node;
            return NULL;
        }
        return node;
    }
    else {
        if (point[dimension] < node -> data[dimension])
            node -> left = removeRec(point, node -> left, depth + 1);
        else 
            node  -> right = removeRec(point, node -> right, depth + 1);   
        return node;
    }
}

void kDTree :: remove(const vector<int> &point){
    this -> root = removeRec(point, this -> root, 0);
}

/*------------------------------------------------------------------*/
/*Eurclid distance*/
double distance(const vector<int> &u, const vector<int> &v){
    // this  two vector must has equally dimension
    assert(u.size() ==  v.size());
    double sum = 0.0;
    double diff;
    for (int i = 0; i < u.size(); i++){
        diff = u[i] - v[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

void kDTree :: nearestNeighbourRec(const vector<int> &target, kDTreeNode *&best, kDTreeNode* node, int depth){
    // R: distance from target to best
    // r: distance from target to current
    // d: distance about dimension alpha from target to current
    double R, r;
    int d;
    /*
    base case:
    node is null : return
    node is leaf and best is null : traversal and choose leaf node for the nearest neighbor currently, return to recursion stack
    */

    if (!node)
        return;
    if (!best && !node -> left && !node ->right){
        best = node;
        return;
    }
    int dimension = depth % k;
    // node can't move: 
    //  - current node equal with tartget choosse it is best and return
    //  - target[dimension] < current[dimension] && leftsubtree is null
    //  - target[dimension] >= current[dimension] && rightsubtree is null
    if (equal(target, node -> data)){
        best = node;
        return;
    }
    // traverse left subtree 
    if (target[dimension] < node-> data[dimension]){
        if (!node -> left) {
            best = node;
            return;
        }
        else 
            nearestNeighbourRec(target, best,node -> left, depth + 1 );
        r = distance(target, node -> data);
        R = distance(target, best -> data);
        d = abs(target[dimension] - node -> data[dimension]);
        if ( r < R) {
            // updatate best 
            best = node;
        }
        // if the plane split cut sphere
        if ( d < R) 
            nearestNeighbourRec(target, best, node -> right, depth + 1);
    }
    // traverse right subtree
    else {
        if (!node -> right){
            best = node;
            return;
        }
        else
            nearestNeighbourRec(target, best,node -> right, depth + 1 );
        r = distance(target, node -> data);
        R = distance(target, best -> data);
        d = abs(target[dimension] - node -> data[dimension]);
        if ( r < R) {
            // updatate best if the distance is shorten than.
            best = node;
        }
        // if the plane split cut sphere
        if ( d < R) 
            nearestNeighbourRec(target, best, node -> left, depth + 1);
    }

}

void kDTree :: nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
    nearestNeighbourRec(target, best, this -> root,0 );
}
int main()
{
    kDTree* T = new kDTree(3);
   int arr[][3] = {{7,5,10}, {9,9,7} ,{10,6,1} ,{5,5,3}, {4,5,1},
                    {1,10,10},{2,8,2} ,{4,1,1}, {1,6,2}, {8,4,5},
                    {3,2,3}, {8,9,6}, {4,7,9}, {4,3,5}, {4,1,1},
                    {10,6,9}, {5,1,7}, {6,5,10}, {7,3,10}, {6,8,2}};
    vector<vector<int>> points;
    int size = sizeof(arr)/sizeof(arr[0]);
    for (int i = 0; i < size; i++){
        vector<int> point;
        for (int j = 0; j < 3; j++)
            point.push_back(arr[i][j]);
        points.push_back(point);
    }
    T -> buildTree(points);
    T -> preorderTraversal();
    cout<<endl;
    T -> inorderTraversal();
    cout<<endl<<"H: "<<T-> height()<<endl;
    cout<<"Leaf: "<<T-> leafCount()<<endl;
    vector<int> target = {10,2,5};
    kDTreeNode* best = NULL;
    T -> nearestNeighbour(target, best);
    T -> visualize();
    cout<<"------------------------------------------------------------------"<<endl;
    best -> print();
    cout<<endl;	
    cout<<"------------------------------------------------------------------"<<endl;

    return 0;
    /*
    x:(2,7) (3,6) (6,12) (9,1) (10,19) (13,15) (17,15) -> (9,1)
    y:(3,6) (2,7) (6,12) -> (2,7)
    y:(13,15) (17,15) (10,19)
    insert (5,6)
    insert (2,2)
    insert (2,8)
    insert (7,3)
    insert (3,5)
    insert (8,1)
    insert (8,7)
    insert (9,6)
    
    // (0,0) (2)
    
    
    
    
        kDTree* T = new kDTree(1);
    int n,d;
    cout<<"number of node && dimension"<<endl;
    cin>>n>>d;
    int x;
    for (int i = 0; i < n; i++){
        cout<<"type node"<<endl;
        vector<int> point;
        for (int j = 0; j < d; j++){
            cin>>x;
            point.push_back(x);
        }
        T -> insert(point);
    }
    T -> inorderTraversal();
    T -> visualize();
    cout<<endl<<"H: "<<T-> height()<<endl;
    cout<<"Leaf: "<<T-> leafCount()<<endl;
    kDTreeNode* best;
    bool flag = true;
    char c;
    while (true){
        cout<<"Type target: "<<endl;
        vector<int> target;
        for (int i = 0; i < d; i++){
            cin >>x;
            target.push_back(x);
        }
        best = nullptr;
        T -> nearestNeighbour(target,best);
        best -> print();
        cout<<"repeat ? y/n"<<endl;
        cin >> c;
        if (c == 'y') flag = true;
        else if (c == 'n') flag = false;
    }
    */
}
