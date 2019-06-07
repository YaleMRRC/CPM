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

