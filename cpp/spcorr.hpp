
/* Implementation of Spearman Rank Correlation Coefficient(SRCC)
Language: c++

Oluwatosin Oluwadare
Department of Computer Science
University of Missouri, Columbia
Columbia
USA

Homepage:https://sites.google.com/site/oluwatosinoe/home
email: 	oeow39@mail.missouri.edu

Requires: pearson correlation coeeficient method (download)
Example:
scc = spearman(vector1, vector2);


Last Update: 2/13/2017

*/

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <iterator> 
#include <vector>
#include <sstream>
#include <random>


using namespace std;


/*
   search for index
   */


void searchList(vector<double> theArray, int sizeOfTheArray, double findFor, vector<int> &index){
	vector<int> foundIndices;

	int j = 0;

	for (int i = 0; i < sizeOfTheArray; i++)
	{
		if (theArray[i] == findFor){
			foundIndices.push_back(i);
			j++;
		}
	}

	if (foundIndices.size()!=0){
		//cout << " Found in index: ";
		for (int i = 0; i < foundIndices.size(); i++){
			// cout << foundIndices[i]+1 << " ";
			index.push_back( foundIndices[i]+1);
		}
	}
	else{
		// cout << " Not found in array";
	}
}

/* Rank the variable vector supplied as input

*/


void Rank(vector<double> vec,vector<double> orig_vect, vector<double> &rank){
	vector<double> R ;
	vector<int> Indices;
	// Declaring new vector and copying element of old vector constructor method, Deep copy
	vector<double> vect2(vec); // vect2 is a sorted list
	int length = vect2.size();
	// assign rank for Sorted list	
	for(int k=0;k<vec.size();k++) {	
		R.push_back(k+1); // starting with 1		
	}	


	// find the unique element in vector
	std::vector<double>::iterator it;
	it = std::unique(vec.begin(), vec.end());
	vec.resize( std::distance(vec.begin(),it) );


	//Break Ties
	for (int k=0;k<vec.size();k++){
		// Search for the index position by value
		Indices.clear();		
		searchList(vect2,length,vec[k],Indices);		
		// Find mean position
		double sum = 0;
		for (int i=0;i<Indices.size();i++){
			sum+=R[Indices[i]-1];
		}				
		double mean_index =   sum / Indices.size();
		//change the rank at ties position
		for(int j=0;j<Indices.size();j++){
			R[Indices[j]-1] = mean_index;
		}		
	}	

	// Search sorted list for index of item on original vector	
	double nPosition;

	for(int k=0; k < orig_vect.size();k++){
		Indices.clear();		
		searchList(vect2,length,orig_vect[k],Indices);
		nPosition = Indices[0]; // just need one ocurrence		
		// Get the respective postion in sorted list then pushback in rank		
		rank.push_back(R[nPosition-1]);		
	}



}

/**
 *	Spearman Correlation
 */
double spearman(vector<double> v1, vector<double> v2){
	comp(v1, v2); // check length
	int n = v1.size();
	vector<double> R1;
	vector<double> R2;
	vector<double> d;
	vector<double> vector1(v1); // original vector v1
	vector<double> vector2(v2); // original vector v2
	// Sort
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());	
	// get the rank on the data
	Rank(v1,vector1, R1); Rank(v2,vector2,R2);
	return pearson(R1,R2);
	//Method 2 : Use the spearman correlation formular( Only if all n ranks are distinct integers )
	for (int k=0;k < n ;k++){		
		double diff = R1[k]- R2[k]; // Difference d where R1.size() = R2.size()
		double sq_diff = pow(diff,2);
		d.push_back(sq_diff );
	}
	// Sum the Squared difference
	double sum = accumulate(d.begin(), d.end(), 0.0);
	int en = n;
	double en3n = (en * en * en) - en;
	double numerator = 6 * sum;
	double corr = 1  - (numerator / en3n);
	return corr;
}

