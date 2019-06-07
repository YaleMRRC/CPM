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
