# Object Oriented C++ CPM

This is object oriented design for connectome based predictive modeling (CPM) by Xilin Shen 2017. Main goals of this
structure are: high interpretability for neuro scientists, capability of fast development, preventing duplicate codes, and eliminating long lines of function parameters which are pretty common in Matlab.
## Class Structure 
This framework has the following base classes:

1. `subject`: each subject has ``` all_edges=K * 268*268``` of connectome where ```K``` is the number of tasks. We also keep the following members: ```num_node```, ```num_task```, and ```id```. 
2. `group`: each consisting of ```N``` subjects. We usually keep group level properties in this class. 
3. `phenotype`: has ```name``` and a list of behavioral measures which is named ```all_behav```.
4. `predictory`: cpm is a predictory model. This is the reason we define two base models : one for predictory models and one for explanatory ones. 
explanatory models don't have k-folds. There is a long history of literature on these models. Generally predictive models produce lower variance and relatively higher bias. 
6. `cpm`: this is the main ridge cpm class with a ```run``` and ```evaluate``` function. 

## K Fold Cross Validation
K fold cross validation is a key module in predictory modules. In order to report full performancer of the CPM model we partition data into `k` folds and use `k-1` fold as training and the `k`th fold as a test.  
To implement K fold cross validation in c++ we first label each element by residual of each index to `k` and then shuffle them. 

## Significant Test
CPM uses a type of feature selection based on t-student test. To this aim, we need to compute correlation between connectivity matrix `x` and behavior vector `y`. We first dmean each of these variables as follows:

```cpp
   double* dmx = dmean(x,n,p1);
   double* dmy = dmean(y,n,p2);
```

where `dmean` subtracts each element from `mean` of the column. If we name `xp` as transpose of matrix `x`, then `coef[j]` which is Pearson correlation is obtained by:

```cpp
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
```
to compute t value and then corrosponding p-value we need to compute beta distribution: ``` beta(z,w) = integral from 0 to 1 of t.^(z-1) .* (1-t).^(w-1) dt. ```. `incbeta` computes incomplete beta distribution.

```cpp
    for(int i=0;i<p1;i++){
                t[i] = coef[i]*pow((double)(n-2.0)/(1.0-pow(coef[i],2)),0.5);
                double a = 0.5,b = (n-2.0)/2.0;
                double c_upper = (double)pow(t[i],2)/(pow(t[i],2)+n);
                double c_lower = (double)(n-2.0)/(pow(t[i],2)+n-2.0);
                if(n>pow(t[i],2)+2){
                        pval[i] =pval[i] = 1.0-this->incbeta(a,b,c_upper);
                }else{
                        pval[i] =pval[i] = this->incbeta(b,a,c_lower);
                }
        }
```
now we can compare `pval` of each edge/feature with predefined threshold on CPM.

## Optimization
In order to make comparisions fast, we can compute inverse of beta distribution and then compae corrosponding `c` values. To do so, we need to first make sure if `t[i]` is in upper of lower part of beta distribution and then save it for later comparision with corrosponding `c` value of threshold. 

```cpp
     double* t = new double[p1];
     double* pval = new double[p1];                                                                
     bool* lower = new bool[p1];                                                                   
     for(int i=0;i<p1;i++){
             t[i] = coef[i]*pow((double)(n-2.0)/(1.0-pow(coef[i],2)),0.5);
             double a = 0.5,b = (n-2.0)/2.0;
             double c_upper = (double)pow(t[i],2)/(pow(t[i],2)+n);
             double c_lower = (double)(n-2.0)/(pow(t[i],2)+n-2.0);
             if(n>pow(t[i],2)+2){
                     pval[i] =c_upper;
                     lower[i]=false;
             }else{
                     pval[i] =c_lower; 
                     lower[i]=true;
             }
     }
```

In the `CPM.cpp` code when we do comparision, we only need to compare `c` values with corrosponding value of the threshold. if `beta(x,a,b)<eta` for `lower` and `1-beta(x,b,a)<eta` for `upper`:
```cpp
      double t1 = xinbta ( b, a, beta_log, this->threshold, ifault );                               
      double t2 = xinbta ( a, b, beta_log, 1.0-this->threshold, ifault ); 
       
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
```

## Prediction
To predict behavior on test fold, we have to fit our training data and appropriate behavioral measures into a line:

```cpp
       double coefficients[order + 1]; // resulting array of coefs
       int result = polyfit(train_sum,
                       ytrain,
                       n-testCount,
                       order,
                       coefficients);
       testInd=0;
       for(int i=0;i<n;i++){
               if(indices[i]==fold){
                       this->predicted[i]=test_sum[testInd]*coefficients[1]+coefficients[0];
                       testInd++;
               }
       }
```

## How to Use
To use ridge cpm (rCPM) you need to load your data in following way:

```cpp
   make -f Makefile
   ./cpm
```

### Build Group
Then you need to build a group of subjects. An example of this function is in ```main.m```: 

```cpp
Group buildGroup(double* phenotype,const group_options opg){                                         
        
        Subject subjects[opg.num_subj];
        for(int i=0;i<opg.num_subj;i++){ 
                phenotype[i] = (rand() % static_cast<int>(20 + 1));                                  
        }
        for(int i=0;i<opg.num_subj;i++){
                double* xi = new double[opg.num_edges];                                              
                for(int j=0;j<opg.num_edges;j++){
                        xi[j] =(rand() % static_cast<int>(4 + 1));                                   
                }       
                subjects[i].setConnectome(xi);                                                       
        }
        return Group(subjects,opg);                                                                  
}
```
and you need to call this with appropriate dataset name and mask matrix over edges. 
If you need to use mask (e.g., abi, nbs) you can add your own one in 
```function this = subject(x,id,dataset,mask)``` in subject.m. Any group level functions can be placed here. 
```C++
    double phenotype[opg.num_subj];
    Group group = buildGroup(phenotype,opg); 
```
CPM needs a number of parameters and we are using Matlab's structures to this aim. ```diagnosis``` is useful when you have group labels.
We initialize it with zero:
```cpp
    cpm_options opc = {};
    group_options opg = {};                                                                      
    opc.threshold = 0.01;                                                                        
    opc.k=3;
    opc.seed=870;
    opc.lambda = 0.0001;                                                                         

    opg.num_task = 0;
    opg.num_node = 268;
    opg.num_edges = 5;//268*268;                                                                 
    opg.num_subj = 10;
```

Finally we need to get an instanciation of cpm and calling ```run()``` and then simply evaluate it. As you see, the constructor in ```rcpm``` only calls ```this = this@predictory(group,options);
```. Other predictive models can include further private members. You can see that function ```run``` in ```predictory``` is ```Abstract```.
This means that all classes that inherits ```predictory``` have to implement it locally.

```cpp
    CPM* c = new CPM(group,phenotype,opc);                                                       
    c->run();
    c->evaluate();
```

