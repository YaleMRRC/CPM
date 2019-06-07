#ifndef RCPM_H
#define RCPM_H
#include "Group.hpp"
#include "predictory.hpp"

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
 * Simulator is a class to simulate the CPM Problem
 *
 */
class RCPM: public predictory{
	public:
		/* constructs
		 * a simulator for numAgents agents. The first numOne of 
		 * these have initial choice 1;
		 * the remainder have initial choice 0. seed is used to i
		 * itialize the random number
		 * generator random().
		 * */
		RCPM(Group group,double* phenotype, cpm_options op);
		/*
		 * runs the simulation for as many rounds as it takesto reach
		 * consensus. The number of communication rounds used is stored in th
		 * output parameter rounds. The consensus value is returned
		 *
		 */
		void run();
};
#endif
