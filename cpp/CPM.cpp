#include "CPM.hpp"
using namespace std;
using namespace Eigen;

/*
 * Here I'm creating a dynamic array of size numAgents
 * I set first numOne agents to be 1 and the others to be zero;
 *
 */
CPM::CPM(Group group,double* phenotype,cpm_options op):predictory(group,phenotype,op){
}

/*
 * Simulation main part. I'm summing all ch values. if it's 1 or 0 then all of them are the same;
 */
void CPM::run(){
	double * X = this->group.getX();
	double * y = this->phenotype;
	int n = this->num_subj;
	int p = this->num_edges;
	double MIN =0.1E-05;
	int order =1;
	int k=this->k;
	double a=0.5,b=(int)((double)(k-1.0)/(k)*n/2-1);//%(n-2)/2.0,xth,fx;
	int ifault;
	//double th = 0.01;//this->threshold;
	//cin>>th;
	//this->beta_inc_values ( n_data, a, b,th, fx );
	double beta_log = lgamma ( a )
		+ lgamma ( b )
		- lgamma ( a + b );
	//cout<<"0.01 -> "<<this->InverseBeta(0.01,0.5,(n-2)/2,0,10)<<" "<<xinbta ( a, b, beta_log,0.01, ifault)<<endl;

	cout<<b<<" "<<a<<" $$$$$ "<<k<<endl;
	double t1 = xinbta ( b, a, beta_log, this->threshold, ifault );
	double t2 = xinbta ( a, b, beta_log, 1.0-this->threshold, ifault );
	

	cout<<"T1:"<<t1<<" T:"<<t2<<endl;//pow(tthresh*(n-2)/(1-tthresh),0.5)<<endl;
	//tthresh = pow(tthresh*(n-2)/(1-tthresh),0.5);

	/*for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			cout<<X[i*p+j]<<" ";
		}
		cout<<endl;
	}

	cout<<"================="<<endl;
	for(int i=0;i<n;i++){
		cout<<y[i]<<" ";
	}
	cout<<endl;
	
	for(int i=0;i<n;i++)
		cout<<indices[i]<<" ";
	cout<<endl;*/
	int* indices = this->kfold(n,this->k);
	for(int fold=0;fold<this->k;fold++){
		cout<<"fold "<<fold<<endl;

		int testCount=0;
		for(int i=0;i<n;i++){
			if(indices[i]==fold){
				testCount++;
			}
		}

		double * xtrain = new double[(n-testCount)*p];
		double * ytrain = new double[n-testCount];

		double * xtest = new double[testCount*p];
		//double * ytest = new double[testCount];

		int trainInd =0;
		int testInd =0;
		for(int i=0;i<n;i++){
			for(int j=0;j<p;j++){
				if(indices[i]==fold){
					xtest[testInd*p+j]=X[i*p+j];
				}else{
					xtrain[trainInd*p+j]=X[i*p+j];
					ytrain[trainInd]=y[i];
				}
			}
			if(indices[i]==fold)
				testInd++;
			else
				trainInd++;
		}

		double * train_sum = new double[n-testCount];
		double * test_sum = new double[testCount];

		for(int i=0;i<n-testCount;i++){
			train_sum[i]=0;
		}

		for(int i=0;i<testCount;i++){
			test_sum[i]=0;
		}


		CORDAT c = this->corr(xtrain,ytrain,n-testCount,p,1);

		for(int i=0;i<n-testCount;i++){
			for(int j=0;j<p;j++){
				if(c.lower[j]){
					if(c.pval[j]<t1 && c.coef[j]>0 ){
						train_sum[i]+=xtrain[i*p+j];
					}
					if(c.pval[j]<t1 && c.coef[j]<0){
						train_sum[i]-=xtrain[i*p+j];
					}
				}else{
					if(c.pval[j]>t2 && c.coef[j]>0 ){
						train_sum[i]+=xtrain[i*p+j];
					}
					if(c.pval[j]>t2 && c.coef[j]<0){
						train_sum[i]-=xtrain[i*p+j];
					}
				}
			}
		}

		for(int i=0;i<testCount;i++){
			for(int j=0;j<p;j++){
				if(c.lower[j]){
					if(c.pval[j]<t1 && c.coef[j]>0 ){
						test_sum[i]+=xtest[i*p+j];
					}
					if(c.pval[j]<t1 && c.coef[j]<0){
						test_sum[i]-=xtest[i*p+j];
					}
				}else{
					if(c.pval[j]>t2 && c.coef[j]>0 ){
						test_sum[i]+=xtest[i*p+j];
					}
					if(c.pval[j]>t2 && c.coef[j]<0){
						test_sum[i]-=xtest[i*p+j];
					}
				}

			}
		}


		double coefficients[order + 1]; // resulting array of coefs
		//		// Perform the polyfit
		int result = polyfit(train_sum,
				ytrain,
				n-testCount,
				order,
				coefficients);

		testInd=0;

		for(int i=0;i<n;i++){
			if(indices[i]==fold){
				this->predicted[i]=test_sum[testInd]*coefficients[1]+coefficients[0];
				if(this->predicted[i]<MIN)
					this->predicted[i] = 0.0;
				testInd++;
			}
		}

		cout<<endl;
		if(result==-1)
			cout<<"polynomial fit before convergance.."<<endl;
		delete[] xtrain;
		delete[] xtest;
		delete[] ytrain;
	}
	
}

int CPM::polyfit(const double* const dependentValues,
		const double* const independentValues,
		unsigned int        countOfElements,
		unsigned int        order,
		double*             coefficients)
{
	// Declarations...
	// ----------------------------------
	enum {maxOrder = 5};

	double B[maxOrder+1] = {0.0f};
	double P[((maxOrder+1) * 2)+1] = {0.0f};
	double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};

	double x, y, powx;

	unsigned int ii, jj, kk;

	// Verify initial conditions....
	// ----------------------------------

	// This method requires that the countOfElements > 
	// (order+1) 
	if (countOfElements <= order)
		return -1;
	// This method has imposed an arbitrary bound of
	// order <= maxOrder.  Increase maxOrder if necessary.
	if (order > maxOrder)
		return -1;
	//cout<<"%%%%"<<endl;

	// Begin Code...
	// ----------------------------------

	// Identify the column vector
	for (ii = 0; ii < countOfElements; ii++)
	{
		x    = dependentValues[ii];
		y    = independentValues[ii];
		powx = 1;

		for (jj = 0; jj < (order + 1); jj++)
		{
			B[jj] = B[jj] + (y * powx);
			powx  = powx * x;
		}
	}

	// Initialize the PowX array
	P[0] = countOfElements;

	// Compute the sum of the Powers of X
	for (ii = 0; ii < countOfElements; ii++)
	{
		x    = dependentValues[ii];
		powx = dependentValues[ii];

		for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
		{
			P[jj] = P[jj] + powx;
			powx  = powx * x;
		}
	}

	// Initialize the reduction matrix
	//
	for (ii = 0; ii < (order + 1); ii++)
	{
		for (jj = 0; jj < (order + 1); jj++)
		{
			A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
		}

		A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
	}

	// Move the Identity matrix portion of the redux matrix
	// to the left side (find the inverse of the left side
	// of the redux matrix
	for (ii = 0; ii < (order + 1); ii++)
	{
		x = A[(ii * (2 * (order + 1))) + ii];
		if (x != 0)
		{
			for (kk = 0; kk < (2 * (order + 1)); kk++)
			{
				A[(ii * (2 * (order + 1))) + kk] = 
					A[(ii * (2 * (order + 1))) + kk] / x;
			}

			for (jj = 0; jj < (order + 1); jj++)
			{
				if ((jj - ii) != 0)
				{
					y = A[(jj * (2 * (order + 1))) + ii];
					for (kk = 0; kk < (2 * (order + 1)); kk++)
					{
						A[(jj * (2 * (order + 1))) + kk] = 
							A[(jj * (2 * (order + 1))) + kk] -
							y * A[(ii * (2 * (order + 1))) + kk];
					}
				}
			}
		}
		else
		{
			// Cannot work with singular matrices
			return -1;
		}
	}

	// Calculate and Identify the coefficients
	for (ii = 0; ii < (order + 1); ii++)
	{
		for (jj = 0; jj < (order + 1); jj++)
		{
			x = 0;
			for (kk = 0; kk < (order + 1); kk++)
			{
				x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
						B[kk]);
			}
			coefficients[ii] = x;
		}
	}

	return 0;
}

