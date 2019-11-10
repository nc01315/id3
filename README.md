# id3
This is a id3 binary classifier (decision tree based) written in C++.
The C++ program split each feature base on the information gain.

The C++ programs (id3.cpp) are developed to creates an ID3 decision tree to classify a set of input training data and then reports the classification performance on a separate set of input testing data. The branch factor of a ID tree is set to 2, in other words , the shape of a ID3 tree is very similar to the shape of a binary search tree , each tree node in ID3 tree contains : a set of 2D data (includes all attributes and class label), a left child pointer and a right child pointers, a double variable splitVal contains an average value of two different adjacent number in the same sorted attribute column, an integer which stores the specific attribute column number the ID3 algorithm is splitting the data set on.


