/*
 * Javid Dadashkarimi
 * jd2392
 *
 */

#include "Subject.hpp"
Subject::Subject(double* x){
	this->x = x;
}

Subject::Subject(){
}

void Subject::setConnectome(double* x){
	this->x = x;
}

double* Subject::getConnectome(){
	return this->x;
}
