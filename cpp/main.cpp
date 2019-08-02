#include <iostream>
#include "tools.hpp"
#include "Group.hpp"
#include "CPM.hpp"
#include <algorithm>
using namespace std;
void run(); // main function
Group buildGroup(double* phenotype,const group_options opg);
/*
 *  aurhor: Javid Dadashkarimi
 *  BioImageSuit Web
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
	opc.k=10;
	opc.seed=870;
	opc.lambda = 0.0001;

	opg.num_task = 0;
	opg.num_node = 268;
	opg.num_edges = 268*267/2;
	opg.num_subj = 50;

	double phenotype[opg.num_subj];
	Group group = buildGroup(phenotype,opg);

	CPM* c = new CPM(group,phenotype,opc); 
	c->run();
	c->evaluate();
}

Group buildGroup(double* phenotype,const group_options opg){
	ifstream inFile1;
	ifstream inFile2;
	inFile1.open("phenotype.txt");
	inFile2.open("connectome.txt");

	Subject subjects[opg.num_subj];
	for(int i=0;i<opg.num_subj;i++){
		//phenfile<< (rand() % static_cast<int>(20 + 1))<<endl;
		string y;
		inFile1>>y;
		//cout<<(rand() % static_cast<int>(20 + 1))<<" ";
		phenotype[i]=stoi(y);//(rand() % static_cast<int>(20 + 1));
	}
	for(int i=0;i<opg.num_subj;i++){
		double* xi = new double[opg.num_edges]; 
		for(int j=0;j<opg.num_edges;j++){
			//myfile<<(rand() % static_cast<int>(4 + 1))<<" ";
			string edge;
			inFile2>>edge;
			//cout<<i<<" "<<j<<" "<<edge<<endl;
			xi[j]=stod(edge);//(rand() % static_cast<int>(4 + 1));
			//cout<<xi[j]<<endl;
		}	
		subjects[i].setConnectome(xi);
	}
	inFile1.close();
	inFile2.close();

	return Group(subjects,opg);
}
