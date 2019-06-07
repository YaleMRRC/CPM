#ifndef Group_H
#define Group_H
#include  "tools.hpp"
#include  "Subject.hpp"
/*
 * @author Javid Dadashkarimi
 * Group is a class to simulate each agent player
 */
struct group_options{
	int num_subj;
	int num_edges;
	int num_task;
	int num_node;
};


class Group{
	private:
		Subject* subjects;
		int group_size;
		int num_node;
		int num_task;
		int num_edges;
		double* X;
	public:
		/*
		 * construct an agent object
		 */
		Group(Subject* subjects, group_options op);
		Group();
		double* getX();
		int getSize();
		int getNumEdges();

};
#endif
