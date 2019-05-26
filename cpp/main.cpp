#include <iostream>
#include "tools.hpp"
#include "Group.hpp"
#include "CPM.hpp"
#include <algorithm>
using namespace std;
void run(); // main function
/*
 *  aurhor: Javid Dadashkarimi
 *  netid: jd2392
 *  Homework: 4
 *  Topic: Consensus
 *  numb of classes: 2 [Simulator, Agent]
 */
int main(int argc, char *argv[]){
	banner();
	run();  
	bye();  
	return 0;
}

void run(){ // main function
	cout<<"cpm "<<endl;
	cpm_options opc = {};
	group_options opg = {};
	opc.threshold = 0.01;
	opc.k=10;
	opc.seed=1000;
	opc.lambda = 0.0001;
	
	opg.num_task = 0;
	opg.num_node = 268;
	opg.num_edges = 10;
	opg.num_subj = 10;
	
	Subject subjects[opg.num_subj];
	double phenotype[opg.num_subj];
	//double* x = new double[num_subj*num_edges]; 

	for(int i=0;i<opg.num_subj;i++){
		phenotype[i] = 30*((double) rand() / (RAND_MAX)) + 1;
	}
	for(int i=0;i<opg.num_subj;i++){
		double* xi = new double[opg.num_edges]; 
		for(int j=0;j<opg.num_edges;j++){
			double a =((double) rand() / (RAND_MAX)) + 1;
			//x[(i-1)*num_edges+j-1] = a;
			xi[j] = a;
		}	
		subjects[i].setConnectome(xi);
	}

	Group group(subjects,opg);
	CPM* c = new CPM(group,phenotype,opc); 
	c->run();
	c->evaluate();
}
