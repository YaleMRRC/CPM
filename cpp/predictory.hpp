#ifndef PREDICTORY_H
#define PREDICTORY_H
#include "Group.hpp"

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
#include <eigen3/unsupported/Eigen/SpecialFunctions>


/*
 * @author: Javid Dadashkarimi
 * Simulator is a class to simulate the Consensus Problem
 *
 */
struct cpm_options{
	float threshold;
	int k;
	int seed;
	float lambda; // value of the lambda, if not provided, cross-validation will be used
	//float v_alpha;
};

struct CORDAT{
	double* coef;
	//double* pval;
	double* pval;
	bool* lower;
};
class predictory{
	protected:
		Group group;
		double* phenotype;
		int num_subj;
		int num_edges;
		float threshold;
		int k;
		int seed;
		float lambda; // value of the lambda, if not provided, cross-validation will be used
		double* predicted;
	public:
		/* constructs
		 * a simulator for numAgents agents. The first numOne of 
		 * these have initial choice 1;
		 * the remainder have initial choice 0. seed is used to i
		 * itialize the random number
		 * generator random().
		 * */
		predictory(Group group,double* phenotype, cpm_options op);
		/*
		 * runs the simulation for as many rounds as it takesto reach
		 * consensus. The number of communication rounds used is stored in th
		 * output parameter rounds. The consensus value is returned
		 *
		 */
		virtual void run()=0;
		CORDAT corr(double* x,double* y,int n,int p1,int p2);
		void evaluate();
		void searchList(vector<double> theArray, int sizeOfTheArray, double findFor, vector<int> &index);
		void Rank(vector<double> vec,vector<double> orig_vect, vector<double> &rank);
		double spearman(vector<double> v1, vector<double> v2);
		double*  dmean(double* x,int n,int p);
		double* vecnorm(double*x,int p1,int p2, int p,int m);
		void tcdf(double*x,double*pval,int n,int p);
		double incbeta(double a, double b, double x);
		int* kfold(int size,int k);
		void shuffleArray(int* array,int size);
		void beta_inc_values ( int &n_data, double &a, double &b, double &x, 
				double &fx );
		double betain ( double x, double p, double q, double beta, int &ifault );
		double r8_max ( double x, double y );
		double xinbta ( double p, double q, double beta, double alpha, int &ifault );
		void timestamp ( );
		double InverseBeta(double p, double alpha, double beta, double A, double B);
};
#endif
