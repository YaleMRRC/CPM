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
	cpm_options opc = {};
	group_options opg = {};
	opc.threshold = 0.01;
	opc.k=3;
	opc.seed=870;
	opc.lambda = 0.0001;
	
	opg.num_task = 0;
	opg.num_node = 268;
	opg.num_edges = 5;//268*268;
	opg.num_subj = 10;
	
	Subject subjects[opg.num_subj];
	double phenotype[opg.num_subj];
	//double* x = new double[num_subj*num_edges]; 

	for(int i=0;i<opg.num_subj;i++){
		phenotype[i] = (rand() % static_cast<int>(20 + 1));
	}

	for(int i=0;i<opg.num_subj;i++){
		double* xi = new double[opg.num_edges]; 
		for(int j=0;j<opg.num_edges;j++){
			xi[j] =(rand() % static_cast<int>(4 + 1));
		}	
		subjects[i].setConnectome(xi);
	}

	Group group(subjects,opg);
	CPM* c = new CPM(group,phenotype,opc); 
	c->run();
	//c->evaluate();
}
