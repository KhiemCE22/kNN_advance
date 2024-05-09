#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */


/*Class kDTreeNode*/



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
/*bonus function for assign operator*/
kDTreeNode* copyTree(kDTreeNode* root){
    if (!root) return NULL;
    kDTreeNode* copy_root = new kDTreeNode(root -> data);
    //clone subtree
    copy_root -> left = copyTree(root -> left);
    copy_root -> right = copyTree(root -> right);
    return copy_root;
}

const kDTree& kDTree :: operator= (const kDTree &other){
    if (this != &other){
        deleteTree(this -> root);
        this -> k = other.k;
        // copy tree
        this -> root = copyTree(other.root);
    }
    return *this;
}

kDTree :: kDTree(const kDTree& other){
    this -> k = other.k;
    this -> root = copyTree(other.root);
}


kDTree :: ~kDTree(){
    deleteTree(this -> root);
}
void kDTree :: inorderRec(kDTreeNode* node) const {
    if (!node) return ;
    inorderRec(node -> left);
    if (node -> left)
        cout<<" ";
    cout<<*node;
    if (node -> right)
        cout<<" ";
    inorderRec(node -> right);
}
void kDTree :: inorderTraversal()const{
    return inorderRec(this -> root);
}
void kDTree :: preorderRec(kDTreeNode* node)const{
    if (!node) return;
    cout<<*node;
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
    cout<<*node;
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


void merge(int left, int middle , int right, vector<vector<int>> &pointList, int dimension, vector<int>& indexOfList){
    // declare leftIndex and rightIndex to save old index
    vector<int> leftIndex, rightIndex;

    vector<vector<int>> leftList(pointList.begin() + left, pointList.begin() + middle + 1);
    vector<vector<int>> rightList(pointList.begin() + middle + 1, pointList.begin() + right + 1);

    if (!indexOfList.empty()){
        leftIndex.assign(indexOfList.begin() + left, indexOfList.begin() + middle + 1);
        rightIndex.assign(indexOfList.begin() + middle + 1, indexOfList.begin() + right + 1);
    }

    int i = 0, j = 0, k = left;
    int n1 = leftList.size();   
    int n2 = rightList.size();

    while (i < n1 && j < n2){
        if (leftList[i][dimension] <= rightList[j][dimension]){
            if (!indexOfList.empty())
                indexOfList[k] = leftIndex[i]; // update  indexOfList
            pointList[k++] = leftList[i++];
        }
        else {
            if (!indexOfList.empty())
                indexOfList[k] = rightIndex[j]; // update  indexOfList
            pointList[k++] = rightList[j++];
        }
    }

    // Tiếp tục ghi đè các phần tử còn lại nếu có
    while (i < n1){
        if (!indexOfList.empty())
            indexOfList[k] = leftIndex[i]; // update  indexOfList
        pointList[k++] = leftList[i++];
    }
    while (j < n2){
        if (!indexOfList.empty())
            indexOfList[k] = rightIndex[j]; // update  indexOfList
        pointList[k++] = rightList[j++];
    }
}


void mergeSort(vector<vector<int>> & pointList, int start, int end, int dimension, vector<int>& indexOfList){
    if (start < end){
        int middle= start + (end - start)/2;
        mergeSort(pointList,start,middle,dimension,indexOfList);
        mergeSort(pointList,middle+1,end,dimension, indexOfList);
        merge(start, middle, end, pointList, dimension, indexOfList);
    }
    
}

void kDTree :: buildTreeRec( vector<vector<int>> &pointList,int start, int end, kDTreeNode* &node, int depth){
    
    vector<int> indexOfList (0);
    if (start > end) return;
    int dimension = depth % k;
    mergeSort(pointList, start, end, dimension, indexOfList);
    int median = start + (end - start)/2;
    

    while(median > start && pointList[median][dimension] == pointList[median-1][dimension]){
        median--;
    }
    vector<int> data = pointList[median];
    node = new kDTreeNode(data);
     // test   
    // cout<<"Create node:"<<endl;
    // cout<<"data: "<<*node<<endl;
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
/*-----------------------------FOR kNN-------------------------------------*/
void kDTree :: buildTreeRec(vector<vector<int>>& pointList, vector<int>& labelList,
                vector<int>&indexOfList, int start, int end, kDTreeNode*& node, int depth){
    if (start > end) return;
    int dimension = depth % k;
    mergeSort(pointList, start, end, dimension, indexOfList);
    int median = start + (end - start)/2;
    
    
    while(median > start && pointList[median][dimension] == pointList[median-1][dimension]){
        median--;
    }
    vector<int> data = pointList[median];
    node = new kDTreeNode(data, labelList[ indexOfList[median] ]);
     // test   
    // cout<<"Create node:"<<endl;
    // cout<<"data: "<<*node<<"label: "<<labelList[ indexOfList[median] ]<<endl;

    buildTreeRec(pointList,labelList, indexOfList,start,median-1,node -> left, depth + 1);
    buildTreeRec(pointList,labelList, indexOfList,median + 1,end,node -> right, depth + 1);   
}
void kDTree :: buildTree(const vector<vector<int>> &pointList, const vector<int>& labelList){
    // destree 
    vector<int> indexOfList;
    for (int i = 0 ; i < pointList.size(); i++)
        indexOfList.push_back(i);
    vector<vector<int>> points = pointList;
    vector<int> labels = labelList;
    if (root)
        deleteTree(this -> root);
    buildTreeRec(points,labels,indexOfList, 0, points.size() - 1, this -> root, 0);

    //test
    // cout<<"Tree with label:"<<endl;
    // cout<<"node: "<<this -> nodeCount()<<endl;
    // cout<<"leaf: "<<this -> leafCount()<<endl;
    // cout<<endl;
    // this -> inorderTraversal();
    // cout<<endl;
    // for (auto i : indexOfList) cout<<i<<" ";
    // cout<<endl;
    // cout<<this -> nodeCount()<<endl;
}
/*-------------------------------------------------------------------------*/
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
    if (!node) return;
    int dimension = depth % this -> k;
    if (target[dimension] < node->data[dimension]){
        nearestNeighbourRec(target, best, node -> left, depth + 1);
        if (!best || distance(target,node->data) < distance(target,best -> data))
            best = node;
        if (abs(node -> data[dimension] - target[dimension]) <  distance(target,best -> data))
            nearestNeighbourRec(target, best, node -> right, depth + 1);
    }
    else {
        nearestNeighbourRec(target, best, node -> right, depth + 1);
        if (!best || distance(target,node->data) < distance(target,best -> data))
            best = node;
        if (abs(node -> data[dimension] - target[dimension]) <  distance(target,best -> data))
            nearestNeighbourRec(target, best, node -> left, depth + 1);
    }
    return;
}
void kDTree :: nearestNeighbour(const vector<int> &target, kDTreeNode *&best){
    best = nullptr;
    nearestNeighbourRec(target, best, this -> root,0 );
}

void last_nullSlot(vector<kDTreeNode*> bestList, int k, int &last, int& nullslot){
    // return the first position equal with NULL 
    nullslot = -1 ;
    for(int i = 0; i < k; i++){
        if (!bestList[i])
            {
                nullslot = i;
                break;
            }
    }
    if (nullslot == -1) last = k - 1; // full
    else if (nullslot == 0) last = -1; // empty
    else last = nullslot - 1; // last = nullslot - 1
    return;
}
void printTestDistance(const vector<int> &target, vector<kDTreeNode *> &bestList){
    for (kDTreeNode* i : bestList){
        if (!i) break;
        cout<<*i<<" d ="<<distance(target, i -> data)<<"|";
    }
    cout<<endl;
}
void updateList(const vector<int> &target, vector<kDTreeNode *> &bestList, kDTreeNode* node, int k){
    double node2tar = distance(target, node -> data);
    for (int i = 0; i < bestList.size(); i++){
        // if NULL Or CLOSER -> insert -> resize -> return 
        if (!bestList[i] || node2tar < distance(target, bestList[i] -> data)){
            bestList.insert(bestList.begin() + i, node);
            if (bestList.size() > k)
                bestList.resize(k);                
            return;
        }
    }
    // if no NULL And No CLOSER And size < k -> push_back
    if (bestList.size() < k){
        bestList.push_back(node);
    }

}



void kDTree :: kNearestNeighbourRec(const vector<int> &target, int k, kDTreeNode* node,vector<kDTreeNode *> &bestList, int depth){

    if (!node) return;

    int dimension = depth % this -> k;
    if (target[dimension] < node -> data[dimension]){
        int last, nullslot;
        kNearestNeighbourRec(target, k, node -> left, bestList, depth + 1);
        last_nullSlot(bestList, k, last, nullslot);
        if ((nullslot != -1) || distance(target, node-> data) < distance(target, bestList[last]->data)){
            updateList(target,bestList,node,k);

        }
        last_nullSlot(bestList, k, last, nullslot);
        if ( (last != -1 && abs(target[dimension] - node->data[dimension]) < distance(target, bestList[last]->data)) || nullslot!= -1)
            kNearestNeighbourRec(target, k, node -> right, bestList, depth + 1);
    }
    else{
        kNearestNeighbourRec(target, k, node -> right, bestList, depth + 1);
        int last, nullslot;
        last_nullSlot(bestList, k, last, nullslot);
        if ((nullslot != -1)|| distance(target, node-> data) < distance(target, bestList[last]->data)){
            updateList(target,bestList,node,k);

        }
        last_nullSlot(bestList, k, last, nullslot);
        if ( (last != -1 && abs(target[dimension] - node->data[dimension]) < distance(target, bestList[last]->data)) || nullslot != -1 )
            kNearestNeighbourRec(target, k, node -> left, bestList, depth + 1);
    }
    return;
}

void kDTree :: kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList){
    int nb_node = this -> nodeCount();

    if (k > nb_node) k = nb_node;
    bestList = vector<kDTreeNode *>(k, nullptr); // create k nullptr fop bestList
    kNearestNeighbourRec(target, k, this -> root, bestList, 0);

}

/*------------------bonus function------------------*/
vector<vector<int>> listlist2vectorvector(list<list<int>> data){
    vector<vector<int>> pointList;
    for (const auto& lst : data) {
        vector<int> vec(lst.begin(), lst.end());
        pointList.push_back(vec);
    }
    return pointList;
}
int major_voting(vector<kDTreeNode*> bestList){
    vector<int> mnist (10,0);
    for (auto i : bestList){
        mnist[i ->label]++;
    }
    int major = mnist[0];
    int predLabel = 0;
    for (int i = 1 ; i < mnist.size(); i++){
        if (major < mnist[i]){
            major = mnist[i];
            predLabel = i;
        }
    }
    return predLabel;
}

/*kNN Class*/
kNN::kNN(int k){
    this -> k = k;
    this -> X_train = nullptr;
    this -> y_train = nullptr;
    this -> Tree = nullptr;
}

void kNN::fit(Dataset& X_train, Dataset& y_train) {
    // No need for dynamic allocation here
    this->X_train = new Dataset(X_train);
    this->y_train = new Dataset(y_train);

    int row, col;
    X_train.getShape(row, col);
    this->Tree = new kDTree(col);   

    vector<vector<int>> pointList = listlist2vectorvector(X_train.data);
    vector<int> labelList;

    for (const list<int>& i : y_train.data) {        
        if (!i.empty()) {
            labelList.push_back(i.front());
        } else {
            throw runtime_error("Empty label encountered");
        }
    }
    Tree->buildTree(pointList, labelList);
    // test

    kDTree* T = new kDTree(this -> k);
    T -> buildTree(pointList);
    // cout<<"Tree without label: "<<endl;
    // cout<<"node: "<<T -> nodeCount()<<"leaf: "<<T ->leafCount()<<endl;
}

Dataset kNN :: predict(Dataset &X_test){
    Dataset y_pred;
    vector<vector<int>> targetList = listlist2vectorvector(X_test.data);
    vector<int> predList;
    for (auto target : targetList){
        vector<kDTreeNode* > bestList (this -> k, nullptr);
        Tree -> kNearestNeighbour(target, this -> k, bestList);
        // major vote
        predList.push_back(major_voting(bestList));
    }
    list<list<int>> labelList;
    for (int i : predList){
        labelList.push_back(list<int>(1,i));
    }
    vector<string> colName (1,"label");
    y_pred.columnName = colName;
    y_pred.data = labelList;
    return y_pred;
}

double kNN::score(const Dataset& y_test, const Dataset& y_pred) {
    list<list<int>> dataTest = y_test.data;
    list<list<int>> dataPred = y_pred.data;
    
    if (dataTest.size() != dataPred.size()) {
        throw runtime_error("Test and prediction data sizes mismatch");
    }

    int sampleSize = dataPred.size();
    if (sampleSize == 0) {
        throw runtime_error("Empty prediction list");
    }

    int count = 0;
    auto itTest = dataTest.begin();
    auto itPred = dataPred.begin();

    for (int i = 0; i < sampleSize; ++i) {
        if (itTest->front() == itPred->front()) {
            count++;
        }
        ++itTest;
        ++itPred;
    }

    return static_cast<double>(count) / sampleSize;
}