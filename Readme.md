## ASSIGNMENT DSA A2 : kNN advance
*This assignment use [kDtree](https://en.wikipedia.org/wiki/K-d_tree) constructure to store multi-dimension data and easy to search k nearest neighbour.*
- **Include files:** 
    - test.cpp : help visualize the kDtree
    - kDtree.hpp kDtree.cpp : implement the class kDtree
    - Dataset.hpp Dataset.o : the header and linker file for class Dataset was build by TA
    - mnist.csv : the mnist dataset, same to DSA A1
- **Use this command to run:**
> g++ -o main main.cpp kDTree.cpp Dataset.o -I . -std=c++11 