#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    int label = -1;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }
    kDTreeNode(vector<int> data, int label, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
        this ->label = label;
    }
    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (int i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
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
    void buildTreeRec(vector<vector<int>>& pointList, vector<int>& labelList, vector<int>& indexOfList, int start, int end, kDTreeNode*& node, int depth);
    // ********************
    kDTreeNode* removeRec(const vector<int> &point, kDTreeNode* & node, int depth);
    void nearestNeighbourRec(const vector<int> &target, kDTreeNode *&best, kDTreeNode* node, int depth);
    void kNearestNeighbourRec(const vector<int> &target, int k, kDTreeNode* node,vector<kDTreeNode *> &bestList, int depth);
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
    void buildTree(const vector<vector<int>> &pointList, const vector<int> &labelList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
    
    // bonus for remove method
    kDTreeNode* getMinNodeRec(kDTreeNode*  node, int split_axis, int depth);
    kDTreeNode* getMinNode(int split_axis);
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    kDTree* Tree;
//    int numClasses;

public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
