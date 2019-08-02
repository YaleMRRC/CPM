#include "predictory.hpp"
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>


using namespace std;

/*
 * Here I'm creating a dynamic array of size numAgents
 * I set first numOne agents to be 1 and the others to be zero;
 *
 */
predictory::predictory(Group group,double* phenotype,cpm_options op){
	this->group = group;
	this->threshold = op.threshold;
	this->k = op.k;
	this->seed = op.seed;
	this->lambda = op.lambda;
	this->phenotype = phenotype;
	this->num_subj = this->group.getSize();
	this->num_edges = this->group.getNumEdges();
	this->predicted = new double[this->num_subj];

}

void predictory::evaluate(){
	cout<<"#########predicted############"<<endl;
	double mse = 0.0;
	double spcorr = 0.0;
	vector<double> y;
	vector<double> yhat;
	for(int i=0;i<this->num_subj;i++){
		y.push_back(this->phenotype[i]);
		yhat.push_back(this->predicted[i]);
		//cout<<this->predicted[i]<<endl;
		mse += pow(this->phenotype[i]-this->predicted[i],2);
	}
	mse /=this->num_subj;
	spcorr = this->spearman(y,yhat);
	cout<<"**** MSE="<<mse<<" *** "<<endl;
	cout<<"**** SPEARMAN="<<spcorr<<" *** "<<endl;
}

CORDAT predictory::corr(double* x,double* y,int n,int p1,int p2){
	string type = "p";
	string rows = "a";
	string tail = "b"; // both
	double* dmx = dmean(x,n,p1);
	double* dmy = dmean(y,n,p2);
	cout<<n<<" "<<p1<<endl;
	double * coef = new double[p1];
	double * xp = new double[n*p1];

	for(int i=0;i<n;i++){
		for(int j=0;j<p1;j++){
			xp[j*n+i]=dmx[p1*i+j];
		}
	}

	for(int j=0;j<p1;j++){
		coef[j] = 0;
		for(int i=0;i<n;i++){
			coef[j]+=xp[j*n+i]*dmy[i];
		}
	}

	double* dx = vecnorm(dmx,n,p1,2,1);
	double* dy = vecnorm(dmy,n,p2,2,1);


	for(int i=0;i<p1;i++){
		coef[i]/=(dx[i]*dy[0]);	
	}

	double* t = new double[p1];
	double* pval = new double[p1];
	bool* lower = new bool[p1];
	for(int i=0;i<p1;i++){
		t[i] = coef[i]*pow((double)(n-2.0)/(1.0-pow(coef[i],2)),0.5);
		//double a = 0.5,b = (n-2.0)/2.0;
		double c_upper = (double)pow(t[i],2)/(pow(t[i],2)+n);
		double c_lower = (double)(n-2.0)/(pow(t[i],2)+n-2.0);
		if(n>pow(t[i],2)+2){
			pval[i] =c_upper;//pval[i] = 1.0-this->incbeta(a,b,c_upper);
			lower[i]=false;
		}else{
			pval[i] =c_lower; //pval[i] = this->incbeta(b,a,c_lower);
			lower[i]=true;
		}
	}

	CORDAT c ={coef,pval,lower};
	return c;
}

	
void predictory::tcdf(double*t,double* pval,int n,int p){
	bool lower =false;
	for(int i=0;i<p;i++){
		if(n<pow(t[i],2)+2)
			lower=true;
	}
	for(int i=0;i<p;i++){
		double a = 0.5;
		double b = n/2.0;
		double c_upper = (double)pow(t[i],2)/(pow(t[i],2)+n);
		double c_lower = (double)n/(pow(t[i],2)+n);
		if(lower)
			pval[i] = this->incbeta(b,a,c_lower);
		else
			pval[i] = 1.0-this->incbeta(a,b,c_upper);
		if(c_lower==0)
			pval[i]=0.0;
	}
	
	return;
}

/*
 *  p indicates norm and m indicates dimension
 */
double* predictory::vecnorm(double*x,int n,int p,int norm,int m){
	if(m==1){
		double* dx = new double[p];
		for(int j=0;j<p;j++){
			dx[j] = 0;
			for(int i=0;i<n;i++){
				dx[j]+=pow(x[i*p+j],2);
			}
			dx[j] = pow(dx[j],1.0/norm);
		}
		return dx;
	}else if(m==2){
		double* dx = new double[n];
		for(int j=0;j<n;j++){
			dx[j] = 0;
			for(int i=0;i<p;i++){
				dx[j]+=pow(x[p*i+j],2);
			}
			dx[j] = pow(dx[j],1.0/norm);
		}
		return dx;
	}else{
		return NULL; 
	}

}

double* predictory::dmean(double* x,int n,int p){
	double* xp = new double[p];
	double* dmx = new double[n*p];

	for(int j=0;j<p;j++){
		xp[j] = 0;
		for(int i=0;i<n;i++){
			xp[j] += x[i*p+j];
			dmx[i*p+j] = x[i*p+j];
		}
		xp[j]/=n;
	}

	for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			dmx[i*p+j] -=xp[j]; 
		}
	}
	delete[] xp;
	return dmx;
}

void predictory::searchList(vector<double> theArray, int sizeOfTheArray, double findFor, vector<int> &index){
	vector<int> foundIndices;

	int j = 0;

	for (int i = 0; i < sizeOfTheArray; i++)
	{
		if (theArray[i] == findFor){
			foundIndices.push_back(i);
			j++;
		}
	}

	if (foundIndices.size()!=0){
		//cout << " Found in index: ";
		for (int i = 0; i < foundIndices.size(); i++){
			// cout << foundIndices[i]+1 << " ";
			index.push_back( foundIndices[i]+1);
		}
	}
	else{
		// cout << " Not found in array";
	}
}

/* Rank the variable vector supplied as input

*/


void predictory::Rank(vector<double> vec,vector<double> orig_vect, vector<double> &rank){
	vector<double> R ;
	vector<int> Indices;
	// Declaring new vector and copying element of old vector constructor method, Deep copy
	vector<double> vect2(vec); // vect2 is a sorted list
	int length = vect2.size();
	// assign rank for Sorted list	
	for(int k=0;k<vec.size();k++) {	
		R.push_back(k+1); // starting with 1		
	}	


	// find the unique element in vector
	std::vector<double>::iterator it;
	it = std::unique(vec.begin(), vec.end());
	vec.resize( std::distance(vec.begin(),it) );


	//Break Ties
	for (int k=0;k<vec.size();k++){
		// Search for the index position by value
		Indices.clear();		
		searchList(vect2,length,vec[k],Indices);		
		// Find mean position
		double sum = 0;
		for (int i=0;i<Indices.size();i++){
			sum+=R[Indices[i]-1];
		}				
		double mean_index =   sum / Indices.size();
		//change the rank at ties position
		for(int j=0;j<Indices.size();j++){
			R[Indices[j]-1] = mean_index;
		}		
	}	

	// Search sorted list for index of item on original vector	
	double nPosition;

	for(int k=0; k < orig_vect.size();k++){
		Indices.clear();		
		searchList(vect2,length,orig_vect[k],Indices);
		nPosition = Indices[0]; // just need one ocurrence		
		// Get the respective postion in sorted list then pushback in rank		
		rank.push_back(R[nPosition-1]);		
	}
}

double predictory::incbeta(double a, double b, double x){

	/*Find the first part before the continued fraction.*/
	const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
	const double beta_ = exp(lbeta_ab);

	//const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;
	int m = 1000;
	double dt = (double)x/m;
	double sum_=1.0e-8;
	for(int i=1;i<m;i++){
		double t = i*dt;
		sum_+= (double)pow(t,a-1)*pow(1-t,b-1)*dt;
	}
	return 1.0/beta_*sum_;
}

/**
 *	Spearman Correlation
 */
double predictory::spearman(vector<double> v1, vector<double> v2){
	//comp(v1, v2); // check length
	int n = v1.size();
	vector<double> R1;
	vector<double> R2;
	vector<double> d;
	vector<double> vector1(v1); // original vector v1
	vector<double> vector2(v2); // original vector v2
	// Sort
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());	
	// get the rank on the data
	Rank(v1,vector1, R1); Rank(v2,vector2,R2);
	//return pearson(R1,R2);
	//Method 2 : Use the spearman correlation formular( Only if all n ranks are distinct integers )
	for (int k=0;k < n ;k++){		
		double diff = R1[k]- R2[k]; // Difference d where R1.size() = R2.size()
		double sq_diff = pow(diff,2);
		d.push_back(sq_diff );
	}
	// Sum the Squared difference
	double sum = accumulate(d.begin(), d.end(), 0.0);
	int en = n;
	double en3n = (en * en * en) - en;
	double numerator = 6 * sum;
	double corr = 1  - (numerator / en3n);
	return corr;
}

int* predictory::kfold(int size,int k){
	int* indices = new int[ size ];

	for (int i = 0; i < size; i++ )
		indices[ i ] = i%k;

	this->shuffleArray( indices, size );

	return indices;
}
void predictory::shuffleArray(int* array,int size) 
{
  int n = size;
  while (n>1) 
  {
    // 0 <= k < n.
    int k = rand()%n;		
    // n is now the last pertinent index;
    n--;					
    // swap array[n] with array[k]
    int temp = array[n];	
    array[n] = array[k];
    array[k] = temp;
  }
}

//****************************************************************************80

void predictory::beta_inc_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
//                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    The incomplete Beta function is also sometimes called the
//    "modified" Beta function, or the "normalized" Beta function
//    or the Beta CDF (cumulative density function.
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//    The function can also be evaluated by using the Statistics package:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 45

  static double a_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.5E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     30.0E+00,
     30.0E+00,
     40.0E+00,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01,
      1.30625,
      1.30625,
      1.30625 };

  static double b_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.0E+00,
      0.5E+00,
      5.0E+00,
      5.0E+00,
     10.0E+00,
      5.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
     11.7562, 
     11.7562, 
     11.7562 };

  static double fx_vec[N_MAX] = {
     0.6376856085851985E-01,
     0.2048327646991335E+00,
     0.1000000000000000E+01,
     0.0000000000000000E+00,
     0.5012562893380045E-02,
     0.5131670194948620E-01,
     0.2928932188134525E+00,
     0.5000000000000000E+00,
     0.2800000000000000E-01,
     0.1040000000000000E+00,
     0.2160000000000000E+00,
     0.3520000000000000E+00,
     0.5000000000000000E+00,
     0.6480000000000000E+00,
     0.7840000000000000E+00,
     0.8960000000000000E+00,
     0.9720000000000000E+00,
     0.4361908850559777E+00,
     0.1516409096347099E+00,
     0.8978271484375000E-01,
     0.1000000000000000E+01,
     0.5000000000000000E+00,
     0.4598773297575791E+00,
     0.2146816102371739E+00,
     0.9507364826957875E+00,
     0.5000000000000000E+00,
     0.8979413687105918E+00,
     0.2241297491808366E+00,
     0.7586405487192086E+00,
     0.7001783247477069E+00,
     0.5131670194948620E-01,
     0.1055728090000841E+00,
     0.1633399734659245E+00,
     0.2254033307585166E+00,
     0.3600000000000000E+00,
     0.4880000000000000E+00,
     0.5904000000000000E+00,
     0.6723200000000000E+00,
     0.2160000000000000E+00,
     0.8370000000000000E-01,
     0.3078000000000000E-01,
     0.1093500000000000E-01,
     0.918884684620518,
     0.21052977489419,
     0.1824130512500673 };

  static double x_vec[N_MAX] = {
     0.01E+00,
     0.10E+00,
     1.00E+00,
     0.00E+00,
     0.01E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     0.50E+00,
     0.90E+00,
     0.50E+00,
     1.00E+00,
     0.50E+00,
     0.80E+00,
     0.60E+00,
     0.80E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.70E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.225609,
     0.0335568,
     0.0295222 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double predictory::betain ( double x, double p, double q, double beta, int &ifault )

//****************************************************************************80
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int &IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
//
{
  double acu = 0.1E-14;
  double ai;
  //double betain;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;

  value = x;
  ifault = 0;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    cerr << "\n";
    cerr << "BETAIN - Fatal error!\n";
    cerr << "  P <= 0.0 or Q <= 0.0\n";
    ifault = 1;
    exit ( 1 );
  }

  if ( x < 0.0 || 1.0 < x )
  {
    cerr << "\n";
    cerr << "BETAIN - Fatal error!\n";
    cerr << "  X < 0.0 or 1 < X\n";
    ifault = 2;
    exit ( 1 );
  }
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }
//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double predictory::r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

void predictory::timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80
double predictory::InverseBeta(double p, double alpha, double beta, double A, double B)
{
	double x = 0;
	double a = 0;
	double b = 1;
	double precision = pow(10, -6); // converge until there is 6 decimal places precision

	while ((b - a) > precision)
	{
		x = (a + b) / 2;
		if (this->incbeta(alpha, beta, x) > p)
		{
			b = x;
		}
		else
		{
			a = x;
		}
	}

	if ((B > 0) && (A > 0))
	{
		x = x * (B - A) + A;
	}
	return x;
}

double predictory::xinbta ( double p, double q, double beta, double alpha, int &ifault )

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    XINBTA computes inverse of the incomplete Beta function.
	//
	//  Discussion:
	//
	//    The accuracy exponent SAE was loosened from -37 to -30, because
	//    the code would not otherwise accept the results of an iteration
	//    with p = 0.3, q = 3.0, alpha = 0.2.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license. 
	//
	//  Modified:
	//
	//    25 September 2014
	//
	//  Author:
	//
	//    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
	//    C++ version by John Burkardt.
	//
	//  Reference:
	//
	//    GW Cran, KJ Martin, GE Thomas,
	//    Remark AS R19 and Algorithm AS 109:
	//    A Remark on Algorithms AS 63: The Incomplete Beta Integral
	//    and AS 64: Inverse of the Incomplete Beta Integeral,
	//    Applied Statistics,
	//    Volume 26, Number 1, 1977, pages 111-114.
	//
	//  Parameters:
	//
	//    Input, double P, Q, the parameters of the incomplete
	//    Beta function.
	//
	//    Input, double BETA, the logarithm of the value of
	//    the complete Beta function.
	//
	//    Input, double ALPHA, the value of the incomplete Beta
	//    function.  0 <= ALPHA <= 1.
	//
	//    Output, int &IFAULT, error flag.
	//    0, no error occurred.
	//    nonzero, an error occurred.
	//
	//    Output, double XINBTA, the argument of the incomplete
	//    Beta function which produces the value ALPHA.
	//
	//  Local Parameters:
	//
	//    Local, double SAE, requests an accuracy of about 10^SAE.
	//
{
	double a;
	double acu;
	double adj;
	double fpu;
	double g;
	double h;
	int iex;
	bool indx;
	double pp;
	double prev;
	double qq;
	double r;
	double s;
	double sae = -30.0;
	double sq;
	double t;
	double tx;
	double value;
	double w;
	double xin;
	double y;
	double yprev;

	fpu = pow ( 10.0, sae );

	ifault = 0;
	value = alpha;
	//
	//  Test for admissibility of parameters.
	//
	if ( p <= 0.0 )
	{
		cerr << "\n";
		cerr << "XINBTA - Fatal error!\n";
		cerr << "  P <= 0.0.\n";
		ifault = 1;
		exit ( 1 );
	}

	if ( q <= 0.0 )
	{
		cerr << "\n";
		cerr << "XINBTA - Fatal error!\n";
		cerr << "  Q <= 0.0.\n";
		ifault = 1;
		exit ( 1 );
	}

	if ( alpha < 0.0 || 1.0 < alpha )
	{
		cerr << "\n";
		cerr << "XINBTA - Fatal error!\n";
		cerr << "  ALPHA not between 0 and 1.\n";
		ifault = 2;
		exit ( 1 );
	}
	//
	//  If the answer is easy to determine, return immediately.
	//
	if ( alpha == 0.0 )
	{
		value = 0.0;
		return value;
	}

	if ( alpha == 1.0 )
	{
		value = 1.0;
		return value;
	}
	//
	//  Change tail if necessary.
	//
	if ( 0.5 < alpha )
	{
		a = 1.0 - alpha;
		pp = q;
		qq = p;
		indx = true;
	}
	else
	{
		a = alpha;
		pp = p;
		qq = q;
		indx = false;
	}
	//
	//  Calculate the initial approximation.
	//
	r = sqrt ( - log ( a * a ) );

	y = r - ( 2.30753 + 0.27061 * r ) 
		/ ( 1.0 + ( 0.99229 + 0.04481 * r ) * r );

	if ( 1.0 < pp && 1.0 < qq )
	{
		r = ( y * y - 3.0 ) / 6.0;
		s = 1.0 / ( pp + pp - 1.0 );
		t = 1.0 / ( qq + qq - 1.0 );
		h = 2.0 / ( s + t );
		w = y * sqrt ( h + r ) / h - ( t - s ) 
			* ( r + 5.0 / 6.0 - 2.0 / ( 3.0 * h ) );
		value = pp / ( pp + qq * exp ( w + w ) );
	}
	else
	{
		r = qq + qq;
		t = 1.0 / ( 9.0 * qq );
		t = r * pow ( 1.0 - t + y * sqrt ( t ), 3 );

		if ( t <= 0.0 )
		{
			value = 1.0 - exp ( ( log ( ( 1.0 - a ) * qq ) + beta ) / qq );
		}
		else
		{
			t = ( 4.0 * pp + r - 2.0 ) / t;

			if ( t <= 1.0 )
			{
				value = exp ( ( log ( a * pp ) + beta ) / pp );
			}
			else
			{
				value = 1.0 - 2.0 / ( t + 1.0 );
			}
		}
	}
	//
	//  Solve for X by a modified Newton-Raphson method,
	//  using the function BETAIN.
	//
	r = 1.0 - pp;
	t = 1.0 - qq;
	yprev = 0.0;
	sq = 1.0;
	prev = 1.0;

	if ( value < 0.0001 )
	{
		value = 0.0001;
	}

	if ( 0.9999 < value )
	{
		value = 0.9999;
	}

	iex = r8_max ( - 5.0 / pp / pp - 1.0 / pow ( a, 0.2 ) - 13.0, sae );

	acu = pow ( 10.0, iex );
	//
	//  Iteration loop.
	//
	for ( ; ; )
	{
		y = betain ( value, pp, qq, beta, ifault );

		if ( ifault != 0 )
		{
			cerr << "\n";
			cerr << "XINBTA - Fatal error!\n";
			cerr << "  BETAIN returned IFAULT = " << ifault << "\n";
			ifault = 1;
			exit ( 1 );
		}

		xin = value;
		y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0 - xin ) );

		if ( y * yprev <= 0.0 )
		{
			prev = r8_max ( sq, fpu );
		}

		g = 1.0;

		for ( ; ; )
		{
			for ( ; ; )
			{
				adj = g * y;
				sq = adj * adj;

				if ( sq < prev )
				{
					tx = value - adj;

					if ( 0.0 <= tx && tx <= 1.0 )
					{
						break;
					}
				}
				g = g / 3.0;
			}
			//
			//  Check whether the current estimate is acceptable.
			//  The change "VALUE = TX" was suggested by Ivan Ukhov.
			//
			if ( prev <= acu || y * y <= acu )
			{
				value = tx;
				if ( indx )
				{
					value = 1.0 - value;
				}
				return value;
			}

			if ( tx != 0.0 && tx != 1.0 )
			{
				break;
			}

			g = g / 3.0;
		}

		if ( tx == value )
		{
			break;
		}

		value = tx;
		yprev = y;
	}

	if ( indx )
	{
		value = 1.0 - value;
	}

	return value;
}
