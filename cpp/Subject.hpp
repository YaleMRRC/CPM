#ifndef Subject_H
#define Subject_H
#include  "tools.hpp"
#include  "Subject.hpp"
/*
 * @author Javid Dadashkarimi
 * Subject is a class to simulate each agent player
 */
class Subject{
	private:
		int group_size;
		int num_node;
		static int num_task;
		static int num_edges;
		int all_edges;
		double* x; 
	public:
		/*
		 * construct an agent object
		 */
		Subject(double* x);
		Subject();
		void setConnectome(double*x);
		double* getConnectome();

};
#endif
