#include "CPM.hpp"
using namespace std;
/*
 * Here I'm creating a dynamic array of size numAgents
 * I set first numOne agents to be 1 and the others to be zero;
 *
 */
CPM::CPM(Group group,double* phenotype,cpm_options op){
	this->group = group;
	this->threshold = op.threshold;
	this->k = op.k;
	this->seed = op.seed;
	this->lambda = op.lambda;
	this->phenotype = phenotype;
}

/*
 * Simulation main part. I'm summing all ch values. if it's 1 or 0 then all of them are the same;
 */
void CPM::run(){
	double ** X = this->group.getX();
}

void CPM::evaluate(){
}

