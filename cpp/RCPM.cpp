#include "RCPM.hpp"
using namespace std;
using namespace Eigen;

/*
 * Here I'm creating a dynamic array of size numAgents
 * I set first numOne agents to be 1 and the others to be zero;
 *
 */
RCPM::RCPM(Group group,double* phenotype,cpm_options op):predictory(group,phenotype,op){
}

/*
 * Simulation main part. I'm summing all ch values. if it's 1 or 0 then all of them are the same;
 */
void RCPM::run(){
	double * X = this->group.getX();
	double * y = this->phenotype;
	CORDAT c = this->corr(X,y,this->num_subj,this->num_edges,1);
	cout<<c.pval[0]<<endl;
}

