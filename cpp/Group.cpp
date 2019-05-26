/*
 * Javid Dadashkarimi
 * jd2392
 *
 */

#include "Group.hpp"
Group::Group(Subject* subjects, group_options op){
	this-> subjects = subjects;
	this-> X =  new double*[this->group_size];
	this->num_node = op.num_node;
	this->num_edges = op.num_edges;
	this->group_size = op.num_subj;
	this->num_task = op.num_task;
	for(int i=0;i<this->group_size;i++){
		X[i] = subjects[i].getConnectome();
	}
}

Group::Group(){
}

double** Group::getX(){
	return this->X;
}
