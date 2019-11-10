//file name : id3.cpp
//Author: Nan Chen

//Object:Develop a software agent in C++ to learn an ID3 decision tree from labeled classification data
//The program should take 3 command-line arguments
// 1) integer: the number of real-valued features in the data set,
// 2) string: input training data filename
// 3) string: input testing data filename
// e.g.
// ./id3 4 training.txt testing.txt
// ./id3 9 trainingCancer.txt testingCancer.txt

//The program should then output the number of testing examples classified correctly by the decision tree

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <numeric>
#include <set>
#include <queue>
using namespace std;

// Attribute sorting
vector<vector<int> > sort_attributes(vector<vector<double> > data)
{
  vector<vector<int> > indices;
  vector<double> *ptr;
  indices.resize(data.size());
  for (int x = 0; x < indices.size(); x++)
	{
    indices[x].resize(data[x].size());
    iota(indices[x].begin(),indices[x].end(),0);
    ptr = &(data[x]);
    sort(indices[x].begin(),indices[x].end(),
	 [&](size_t i, size_t j){ return (*ptr)[i] < (*ptr)[j]; });
  }
  return indices;
}

struct split_node
{
  double IG; 	// information gain
  double L2NormSquaredDistance; // L2-Norm-Squared Distance
  int Row;    // Spliting row info
  int Col;    // Spliting Attribute_i info
  double splitDoubleValue; 	// split

	split_node() {} // default constructor
	bool operator < (const split_node& Rhs) const // Max Heap
	{
		if (IG < Rhs.IG) return true;
		else if (IG == Rhs.IG)            // IG tie handling
    {
        if (L2NormSquaredDistance > Rhs.L2NormSquaredDistance) return true;
        else return false;
    }
    else return false;
	}
}; // end of struct node

// The node class will represent each node in the ID3 tree
struct Node
{
  vector<vector<double> > dataSet; // data set contains all atrributes and class labels
  vector<vector<double> > dataSetL; // upper subset of data set , x <= sp.
  vector<vector<double> > dataSetR; // lower subset of data set , x > sp

  double IG;                       // max information gain
  double Dist;            // L2 square distance
  int rowNum;          // split row position (attribute)
  int colNum;        // split col position
  Node*  leftNodePtr;		  // left Child of a ID3 tree node
  Node*  rightNodePtr;		// right Child of a ID3 tree node
  double splitVal;        // floating splitting number
  vector<int> NodeLabel;              // label for data set
  bool terminalBoolState;   // indidates the node is a termal node or not
  int label;            // 0 , 1, 2 ,3 ...n  , integer class label of node

  Node(){               // default constructor
	IG = 0.0;
	Dist = 0.0;
    leftNodePtr = NULL;
    rightNodePtr = NULL;
    splitVal = 0.0;
    terminalBoolState = false;
    label = -1;
  }

  void FindSlit ()
  {
        priority_queue< split_node > mypq;  // split_node with highest IG on the top of mypq
        double currentVal;
        int currentPos;
        double nextVal;
        double avergeSplit;
        set<int> uniLabelSet; // unique int label in the bag
        double I;              // total Information I(C1,C2...Cn)
        double I_LessEq;
        double I_Greater;
        double E;               // Expectation
        double IG;              // Information Gain (Attribute_i)
        double sum = 0.0;              //
        double sum2 = 0.0;              //
        double sum3 = 0.0;              //
        double P_lessThanEqual;     // P(x <= sp.)
        double P_greaterThan;       // P(x > sp.)
        int splitPos;        // indicates the position of the atrributes
        int ii;
        split_node split_node1;   // highest IG
        split_node split_node2;   //  highest IG
        vector<vector<int> > indices;
        vector< vector< double > > SortedDataSet (dataSet.size(), std::vector<double> (dataSet[0].size(), 0 ) );

    	indices = sort_attributes(dataSet);

        for (int k = 0; k < 1; k++)
        {
    	    //cout << "Sorting by " << k << "th row..." << endl;
    	    
			for (int i = 0; i < dataSet[0].size(); i++)   // col = 1
            {
			  for (int j = 0; j < dataSet.size(); j++)  // row = 5
              {
              SortedDataSet [j][i] = dataSet[j][indices[k][i]] ;
//              cout << SortedDataSet [j][i]<< " ";
    	      }//for
    	     // 	cout << endl;
    	    } //for
    	  }//for

        ii =  dataSet.size()-1;  // maxRowNum -1
  //      cout << "Last Row pos = " << ii <<endl;

        for (int j = 0; j < dataSet[0].size(); j++)  // Col
        {
              uniLabelSet.insert(SortedDataSet [ii][j]);   // class label set
        } //for

    //    cout << "uniLabelSet size = " <<uniLabelSet.size()<<endl;
    //    cout << "col = " <<dataSet[0].size()<<endl;
        vector <int> N( uniLabelSet.size() ,0); // count of data set with unique class label
        vector <int> N_LessEq( uniLabelSet.size() ,0); // count of data set with unique class label , given x <= sp.
        vector <int> N_Greater( uniLabelSet.size() ,0); // count of data set with unique class label , given x > sp.
        vector <double> P( uniLabelSet.size() ,0.0); // prior probabilies for each class label
        vector <double> P_LessEq( uniLabelSet.size() ,0.0); //
        vector <double> P_Greater( uniLabelSet.size() ,0.0); //
        vector <double> DistanceVector ; // Degree of best prior probabilies preservation
        double Distance = 0.0 ; // Degree of best prior probabilies preservation
        set<int>::iterator iter;
        int myDisIndex;

        //i =  dataSet.size()-1;  // maxColNum -1
        for (int i = 0 ; i < N.size();i++) N[i]=0; // reset a vector to all zeros

        for (int j = 0; j < dataSet[0].size(); j++)  // col = 1
        {
            iter = uniLabelSet.find(SortedDataSet [ii][j]);
            if ( iter!= uniLabelSet.end() )
            {
              myDisIndex = distance(uniLabelSet.begin(), iter);
              N[myDisIndex]++;
            } //if
        } //for


        // element in set
        //for (auto itr = uniLabelSet.begin();itr != uniLabelSet.end(); itr++) cout << "In Set: " <<  *itr <<endl;

        // N
        //for (auto itr = N.begin();itr != N.end(); itr++) cout << "N : " <<  *itr <<endl;


        // P
        sum = 0.0 ;
        for (int k = 0; k < uniLabelSet.size(); k++)  // uniLabelSet.size() == 1
        {
            P[k] = (double)N[k]/(double)dataSet[0].size();
            if (P[k] == 0.0 ) continue;
            else
            sum = sum + P[k]* (-1*log2(P[k]));
        } //for

        I = sum ; // total information I (C1,C2,...Cn)
      //  cout << "I = " << I <<endl;

		if (dataSet[0].size() > 1) // DS.size() >= 2
		{
			// we can split data set 
			//i = 0 ; //first row

		//	cout << "the size of sub DS = " << dataSet.size() << endl;
		//	cout << "the size of sub DS = " << dataSet[0].size() << endl;

			for (int iii = 0; iii < dataSet.size() - 1; iii++)   // all row , except for the last row of class label
			{
			//	cout << "sorting " << iii+1 <<" th attribute" << endl;

				for (int k = iii; k < iii + 1; k++)  // sort n th atrribute
				{
					for (int i = 0; i < dataSet[0].size(); i++)   // col = 2
					{
						for (int j = 0; j < dataSet.size(); j++)  // row = 5
						{
							SortedDataSet[j][i] = dataSet[j][indices[k][i]];
				//			cout << SortedDataSet[j][i] << " ? ";
						}
					//	cout << endl;
					} //for
				}//for

				//cout << "Done sorting " << iii + 1 << " th attribute" << endl;

				for (int j = 0; j < dataSet[0].size() - 1; j++)  // col
				{
					currentVal = SortedDataSet[iii][j];
					nextVal = SortedDataSet[iii][j + 1];     //SortedDataSet [j+1][i];
					if (currentVal == nextVal && j+1 < dataSet[0].size() - 1)   // crash in [5][2] case
						continue;
					else if (iii != (dataSet.size() - 1) && currentVal == nextVal && j + 1 >= dataSet[0].size() - 1)
					{
						break;
					}
					else if (iii == (dataSet.size() - 1) && currentVal == nextVal && j + 1 >= dataSet[0].size() - 1)
					{
						// we no longer can split data set , terminal node , assign class label
						terminalBoolState = true;
						label = dataSet[dataSet.size() - 1][0];   // pick the 1st data's label as node's label 
						leftNodePtr = NULL;   // Left Child Ptr NULL
						rightNodePtr = NULL;    // Right Child Ptr NULL
						return;
					}
					else
					{
						splitPos = j;   // store aP_LessEqrox. slit row position
						avergeSplit = (currentVal + nextVal) / 2.0;
						//cout << "currentVal = " << currentVal << endl;
						//cout << "nextVal = " << nextVal << endl;
						//cout << "avergeSplit = " << avergeSplit << endl;
						P_lessThanEqual = (double)(j + 1) / (double)dataSet[0].size();
						P_greaterThan = 1.0 - P_lessThanEqual;

						// Left part in a column , x <= sp.
						for (int i = 0; i < N_LessEq.size(); i++) N_LessEq[i] = 0; // reset a vector to all zeros

						for (int j = 0; j <= splitPos; j++)  // row
						{
							iter = uniLabelSet.find(SortedDataSet[ii][j]);
							if (iter != uniLabelSet.end())
							{
								myDisIndex = distance(uniLabelSet.begin(), iter);
								N_LessEq[myDisIndex]++;
							} //if
						} //for

						  // P_LessEq
						sum2 = 0.0;
						for (int k = 0; k < uniLabelSet.size(); k++)  //
						{
							P_LessEq[k] = (double)N_LessEq[k] / (double)(splitPos + 1);
							if (P_LessEq[k] == 0.0) continue;
							else sum2 = sum2 + P_LessEq[k] * (-1 * log2(P_LessEq[k]));
						} //for
						I_LessEq = sum2;

						// Right part in a column ,  x > sp.
						for (int i = 0; i < N_Greater.size(); i++) N_Greater[i] = 0; // reset a vector to all zeros

						for (int j = splitPos + 1; j < dataSet[0].size(); j++)  // row
						{
							iter = uniLabelSet.find(SortedDataSet[ii][j]);
							if (iter != uniLabelSet.end())
							{
								myDisIndex = distance(uniLabelSet.begin(), iter);
								N_Greater[myDisIndex]++;
							} //if
						} //for

						  // P_Greater
						sum3 = 0.0;
						for (int k = 0; k < uniLabelSet.size(); k++)  //
						{
							P_Greater[k] = (double)N_Greater[k] / (double)(SortedDataSet[0].size() - splitPos - 1);
							if (P_Greater[k] == 0.0) continue;
							else sum3 = sum3 + P_Greater[k] * (-1 * log2(P_Greater[k]));
						} //for

						I_Greater = sum3;

						E = P_lessThanEqual * (I_LessEq)+P_greaterThan * (I_Greater);

						for (int k = 0; k < uniLabelSet.size(); k++)  // L2 Distance
						{
							Distance = Distance + pow(P[k] - P_LessEq[k], 2.0) + pow(P[k] - P_Greater[k], 2.0);
						} //for

						DistanceVector.push_back(Distance);  // store all L2 Distance in a vector ?

						//calculate information gain = I(C1...Cn) - E(Attribute_i)
						IG = I - E;
						//cout << "P_lessThanEqual = " << P_lessThanEqual << endl;
						//cout << "P_greaterThan = " << P_greaterThan << endl;
						//cout << "I_LessEq = " << I_LessEq << endl;
						//cout << "I_Greater = " << I_Greater << endl;
						//cout << "E = " << E << endl;
						//cout << "Distance = " << Distance << endl;
						//cout << "Number of Splits = " << DistanceVector.size() << endl;
						//cout << "IG = " << IG << endl << endl;

						split_node split_nodeX;
						split_nodeX.IG = IG;
						split_nodeX.L2NormSquaredDistance = Distance;
						split_nodeX.Row = iii;
						split_nodeX.Col = j;
						split_nodeX.splitDoubleValue = avergeSplit;
						mypq.push(split_nodeX);
					} //END else
				} //for
			} //for

			for (auto it = uniLabelSet.begin(); it != uniLabelSet.end(); ++it)
				NodeLabel.push_back(*it);

			if (mypq.size() >= 2)
			{
				split_node1 = mypq.top();			// get the split with highest IG // crash !
				mypq.pop();				// pop the highest IG node

				split_node2 = mypq.top();			// get the split with highest IG
				mypq.pop();				// pop the highest IG node

			//split_node1 and split_node2 may have the same IG and L2 distance
				if (split_node1.IG > split_node2.IG)
				{
					IG = split_node1.IG;
					Dist = split_node1.L2NormSquaredDistance;
					rowNum = split_node1.Row;       // split row position
					colNum = split_node1.Col;        // split col position
					splitVal = split_node1.splitDoubleValue;  // floating splitting number
				}

				else if (split_node1.IG < split_node2.IG)
				{
					IG = split_node2.IG;
					Dist = split_node2.L2NormSquaredDistance;
					rowNum = split_node2.Row;       // split row position
					colNum = split_node2.Col;        // split col position
					splitVal = split_node2.splitDoubleValue;  // floating splitting number
				}
				else // split_node1.IG == split_node2.IG , tie
				{
					if (split_node1.L2NormSquaredDistance < split_node2.L2NormSquaredDistance)
					{
						IG = split_node1.IG;
						Dist = split_node1.L2NormSquaredDistance;
						rowNum = split_node1.Row;       // split row position
						colNum = split_node1.Col;        // split col position
						splitVal = split_node1.splitDoubleValue;  // floating splitting number
					}
					else if (split_node1.L2NormSquaredDistance > split_node2.L2NormSquaredDistance)
					{
						IG = split_node2.IG;
						Dist = split_node2.L2NormSquaredDistance;
						rowNum = split_node2.Row;       // split row position
						colNum = split_node2.Col;        // split col position
						splitVal = split_node2.splitDoubleValue;  // floating splitting number
					}
					else // split_node1.L2NormSquaredDistance == split_node2.L2NormSquaredDistance , tie
					{
						if (split_node1.Row < split_node2.Row)// pick the split comes first
						{
							IG = split_node1.IG;
							Dist = split_node1.L2NormSquaredDistance;
							rowNum = split_node1.Row;       // split row position
							colNum = split_node1.Col;        // split col position
							splitVal = split_node1.splitDoubleValue;  // floating splitting number
						}
						else if (split_node1.Row > split_node2.Row)// pick the split comes first
						{
							IG = split_node2.IG;
							Dist = split_node2.L2NormSquaredDistance;
							rowNum = split_node2.Row;       // split row position
							colNum = split_node2.Col;        // split col position
							splitVal = split_node2.splitDoubleValue;  // floating splitting number
						}
						else // split_node2 and split_node1 on the same Row
						{
							if (split_node1.Col < split_node2.Col)
							{
								IG = split_node1.IG;
								Dist = split_node1.L2NormSquaredDistance;
								rowNum = split_node1.Row;       // split row position
								colNum = split_node1.Col;        // split col position
								splitVal = split_node1.splitDoubleValue;  // floating splitting number
							}
							else
							{
								IG = split_node2.IG;
								Dist = split_node2.L2NormSquaredDistance;
								rowNum = split_node2.Row;       // split row position
								colNum = split_node2.Col;        // split col position
								splitVal = split_node2.splitDoubleValue;  // floating splitting number
							} // else
						} //else
					} //else
				} //else

			}// if
			else // mypq.size < 2
			{
				terminalBoolState = true;
				label = dataSet[dataSet.size() - 1][0];
				leftNodePtr = NULL;   // Left Child Ptr NULL
				rightNodePtr = NULL;    // Right Child Ptr NULL
				return;
			}
		
			//cout << "Highest IG = " << IG << endl;
			//cout << "L2NormSquaredDistance = " << Dist << endl;
			//cout << "rowNum = " << rowNum << endl;
			//cout << "colNum = " << colNum << endl;
			//cout << "splitVal = " << splitVal << endl;

			if (DistanceVector.size() == 0)  // terminal node
			{
				// assign node with a class label
				terminalBoolState = true;
			}
			else
			{
				terminalBoolState = false;
			}

			if (terminalBoolState == true)
			{
				// assign class label with majority vote
				// if two class labels have a tie , choose the 1st class label
				auto itr = max_element(begin(N), end(N)); // c++11 , max_element() find the 1st max item and return its iterator
				int index = distance(N.begin(), itr);
				label = NodeLabel[index];
				leftNodePtr = NULL;   // Left Child Ptr NULL
				rightNodePtr = NULL;    // Right Child Ptr NULL
			}

			int i = rowNum;   // the atrribute has the Highest IG

			for (int k = i; k < i + 1; k++)	// re-sorting the data set base on the attribute with Highest IG 
			{
				//cout << "Sorting by " << k << "th row..." << endl;

				for (int i = 0; i < dataSet[0].size(); i++)   // col = 2
				{
					for (int j = 0; j < dataSet.size(); j++)  // row = 5
					{
						SortedDataSet[j][i] = dataSet[j][indices[k][i]];
					//	cout << SortedDataSet[j][i] << " ";
					}
				//	cout << endl;
				} //for
			}//for

			dataSetL.resize(dataSet.size());  // determine the row size = 5 of dataSetL
			for (int i = 0; i < dataSet.size(); ++i)
			{
				dataSetL[i].resize(colNum + 1);// determine the col size = 3 of dataSetL
			}

			dataSetR.resize(dataSet.size());  // determine the row size = 5 of dataSetR
			for (int i = 0; i < dataSet.size(); ++i)
			{
				dataSetR[i].resize(dataSet[0].size() - colNum - 1);// determine the col size = 97 of dataSetR
			}

			//cout << "dataSet.size() = " << dataSet.size() << endl;  //5

			// seperate the SortedDataSet into left part , x <= sp.
			for (int i = 0; i < dataSet.size(); i++)  // row
			{
				for (int j = 0; j <= colNum; j++)   // col
				{
					dataSetL[i][j] = SortedDataSet[i][j]; // store the left half
				}
			}

			// seperate the SortedDataSet into right part , x > sp.
			for (int i = 0; i < dataSet.size(); i++)  // row
			{
				for (int j = colNum + 1; j < dataSet[0].size(); j++)   // col
				{
					//cout << "Begin" << endl;
					//cout << "i = " << i << endl;
					//cout << "j = " << j << endl;
					dataSetR[i][j - colNum - 1] = SortedDataSet[i][j]; // store the lower half
				 //cout << "done" << endl << endl;

				}
			}
		}

		else
		{
			// we no longer can split data set , terminal node , assign class label
			terminalBoolState = true;
			label = dataSet[dataSet.size()-1][0];
			leftNodePtr = NULL;   // Left Child Ptr NULL
			rightNodePtr = NULL;    // Right Child Ptr NULL
			return;
		}
 
  } // end of FindSlit ()

}; // end of struct

Node* ID3BinarySplit(vector<vector<double> > data)
{
	//cout << "Input data Row size " << data.size() << endl;
	//cout << "Input data Col size " << data[0].size() << endl << endl;

	Node *newnode = new Node;
	newnode->dataSet = data;
	Node *current = newnode;
	
	//cout << "--------Before Splitting --------" << endl;
	newnode->FindSlit();
	//cout << "--------After Splitting --------" << endl;


	//cout << "Left Splitted data Row size " << newnode->dataSetL.size()<< endl;
	//cout << "Right Splitted Col size " << newnode->dataSetL[0].size() << endl << endl;
	//cout << "Left Splitted data Row size " << newnode->dataSetR.size() << endl;
	//cout << "Right Splitted Col size " << newnode->dataSetR[0].size() << endl << endl;

	if (newnode->terminalBoolState == false) // not termial node 
	{
		newnode->leftNodePtr = ID3BinarySplit(newnode->dataSetL);
		newnode->rightNodePtr = ID3BinarySplit(newnode->dataSetR);
	}

	else
		return current;  // base-case: terminal node , done with splitting
	return current;
}

static int correctlyClassifiedDataCounter;  //????

void traverseID3Tree(vector<double>  dataSet, Node *NodePtr)
{
	int TreeNodeLabel;
	int attributeNum;
	double temp;

	if (NodePtr->terminalBoolState == true) // this node is terminal node 
	{
		TreeNodeLabel = NodePtr->label;
		if (TreeNodeLabel == dataSet[dataSet.size() - 1]) // two labels match
		{
			correctlyClassifiedDataCounter++;  // define to Global ?
			return;
		}
		else
		{
			return; //continue;
		}
	} //if

	else  // this node is an internal node
	{
		attributeNum = NodePtr->rowNum;		// retrieve attributeNum
		temp = dataSet[attributeNum];
		if (temp <= NodePtr->splitVal)	//recursive call myself with LC: NodePtr->leftNodePtr
			traverseID3Tree(dataSet, NodePtr->leftNodePtr);
		else // temp > head->splitVal ,recursive call myself with RC: NodePtr->rightNodePtr
			traverseID3Tree(dataSet, NodePtr->rightNodePtr);
	} //else

	return;
}  // traverseID3Tree ()

int ID3Classification(vector<vector<double> > dataSet, Node *head)
{
	int attributeNum;
	int TreeNodeLabel;
	int temp;
	vector<double> tempVector;
	vector<vector<int> > indices;
	vector< vector< double > > SortedDataSet(dataSet.size(), std::vector<double>(dataSet[0].size(), 0));
	indices = sort_attributes(dataSet);

	for (int k = 0; k < 1; k++)
	{
		//cout << "Sorting by " << k << "th row..." << endl;
		for (int i = 0; i < dataSet[0].size(); i++)   // col = 4
		{
			for (int j = 0; j < dataSet.size(); j++)  // row = 5
			{
				SortedDataSet[j][i] = dataSet[j][indices[k][i]];
				//cout << SortedDataSet[j][i] << " ";
			}//for
			//cout << endl;
		} //for
	}//for

	//... ...
	tempVector.resize(SortedDataSet.size() );
	for (int i = 0; i < SortedDataSet[0].size(); i++)   // col = 4
	{
		tempVector.clear();
		for (int j = 0; j < SortedDataSet.size(); j++)  // row = 5
		{
			tempVector.push_back(SortedDataSet[j][i]);  ;
		}
		traverseID3Tree(tempVector, head);
	}

	return correctlyClassifiedDataCounter;
}  //ID3Classification



int main(int argc, char* argv[])
{
	if (argc != 4)  // must have 3 arguments
	    {
	      cerr << endl;
	      cerr << "Usage: " << argv[0] << " [integer of total number attributes ] [training data file] [testing data file]" << endl;
	      cerr << endl;
	      return 1;
	    }

	  int numAttributes = 4;  // total number attributes of data set
	  string trainingData = argv[2]; // trainingData file name; 
	  string testingData = argv[3]; // testingData file name; 
	  
		//string trainingData = "training.txt"; // trainingData file name
		//string trainingData = "trainingCancer.txt"; // trainingData file name

		//string testingData = "testing.txt";  // testingData file name
		//string testingData = "testingCancer.txt"; // trainingData file name


		ifstream infileTraining;  // training 
		ifstream infileTesting;		// testing 

		infileTraining.open(trainingData);
		infileTesting.open(testingData);

		vector<vector<double> > data;
		vector<vector<double> > data2;

		vector<vector<int> > indices;
		vector<vector<int> > indices2;

		string line;
		string line2;

		double value;
		double value2;

	//Process training file 
		getline(infileTraining,line);
	  stringstream parsed(line);

	  // Prep vectors...
	  while (!parsed.eof()) {
	    parsed >> value;
	    data.push_back(vector<double>());
	  }

	  while (! infileTraining.eof()) {
	    stringstream parsed(line);
	    for (int i = 0; i < data.size(); i++){
	      parsed >> value;
	      data[i].push_back(value);
	    }
	    getline(infileTraining,line);
	  }

	  //Process testing file 
	  getline(infileTesting, line2);
	  stringstream parsed2(line2);

	  // Prep vectors...
	  while (!parsed2.eof()) {
		  parsed2 >> value2;
		  data2.push_back(vector<double>());
	  }

	  while (!infileTesting.eof()) {
		  stringstream parsed2(line2);
		  for (int i = 0; i < data.size(); i++) {
			  parsed2 >> value2;
			  data2[i].push_back(value2);
		  }
		  getline(infileTesting, line2);
	  }
	  
	  Node *head = ID3BinarySplit(data);		// recursively building ID3 tree
	  int correctlyClassifiedDataCounter = ID3Classification(data2, head);
	 // cout << "Num of Data correcly Classified = " << correctlyClassifiedDataCounter << endl;
	cout << correctlyClassifiedDataCounter << endl;

	infileTraining.close();
    infileTesting.close();
  return 0;
}
