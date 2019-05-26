#ifndef CPM_H
#define CPM_H
#include "Group.hpp"
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
	float v_alpha;
};

class CPM{
	private:
        Group group;
	double* phenotype;
        float v_alpha; /* value of the alpha parameter in elastic
        //                                   net, default is 1e-6 which makes the
        //                                   regression method to be ridge
        //                                   regression, v_alpha=1 makes it lasso. */
        int num_subj;
        int num_node;
        int num_task;
        int num_edges;
	float threshold;
	int k;
	int seed;
        float lambda; // value of the lambda, if not provided, cross-validation will be used

	public:
		/* constructs
		 * a simulator for numAgents agents. The first numOne of 
		 * these have initial choice 1;
		 * the remainder have initial choice 0. seed is used to i
		 * itialize the random number
		 * generator random().
		 * */
		CPM(Group group,double* phenotype, cpm_options op);
		/*
		 * runs the simulation for as many rounds as it takesto reach
		 * consensus. The number of communication rounds used is stored in th
		 * output parameter rounds. The consensus value is returned
		 *
		 */
		void run();
		void evaluate();
};
#endif
