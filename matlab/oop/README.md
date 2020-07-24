# Object Oriented ridge CPM

This is object oriented design for connectome based predictive modeling (CPM) by Xilin Shen 2017. Main goals of this structure are
to be easy to understand for neuro scientists, capability of fast development, preventing duplicate codes, and eliminating long lines of function parameters which are pretty common in Matlab.
## Class Structure 
This framework has the following base classes:

1. `subject`: each subject has ``` all_edges=K * 268*268``` of connectome where ```K``` is the number of tasks. We also keep the following members: ```num_node```, ```num_task```, and ```id```. 
2. `group`: each consisting of ```N``` subjects. We usually keep group level properties in this class. 
3. `phenotype`: has ```name``` and a list of behavioral measures which is named ```all_behav```.
4. `predictory`: cpm is a predictory model. This is the reason we define two base models : one for predictory models and one for explanatory ones. 
explanatory models don't have k-folds. There is a long history of literature on these models. Generally predictive models produce lower variance and relatively higher bias. 
5. `explanatory`: this is designed for explanatory models without k-fold. In future versions of this framework, a number of models including but not limited to cca and manova
will inherit from this class. 
6. `rcpm`: this is the main ridge cpm class with a ```run``` and ```evaluate``` function. 

## How to Use
To use ridge cpm (rCPM) you need to load your data in following way:

```Matlab
    dataset = "hcp.50"; % LDA on UCLA + ages on 3 bins for HCP
    x = load('../data.50/HCP900_rest_n50.mat');
    y = load('../data.50/HCP900_PMAT24_A_CR_n50.mat');
    x= x.HCP900_rest_n50; % makes sure x is 268*268*N
    y=y.HCP900_PMAT24_A_CR_n50; % makes sure y is N*1
    N = size(x,3);
```
### Build Group
Then you need to build a group of subjects. An example of this function is in ```main.m```: 

```Matlab
function g = buildGroup(x,dataset,mask)
    N =size(x,3);
    subjects(1,N) = subject(N);
    for i=1:N
        subjects(i) = subject(x(:,:,i,:),i,dataset,mask);
    end
    g = group(subjects);
end
```
and you need to call this with appropriate dataset name and mask matrix over edges. 
If you need to use mask (e.g., abi, nbs) you can add your own one in 
```function this = subject(x,id,dataset,mask)``` in subject.m. Any group level functions can be placed here. 
```Matlab
    g = buildGroup(x,dataset,'none'); % mask=false, Bins
```
CPM needs a number of parameters and we are using Matlab's structures to this aim. ```diagnosis``` is useful when you have group labels.
We initialize it with zero:
```Matlab
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = 2;
    options.phenotype = phenotype('behav',y);
    options.diagnosis = zeros(N,1);
```

Finally we need to get an instanciation of cpm and calling ```run()``` and then simply evaluate it. As you see, the constructor in ```rcpm``` only calls ```this = this@predictory(group,options);
```. Other predictive models can include further private members. You can see that function ```run``` in ```predictory``` is ```Abstract```.
This means that all classes that inherits ```predictory``` have to implement it locally.
```Matlab
    m = cpm(g,options);
    m.run();
    m.evaluate();
```
## Example (manova)
As an example lets take a look at multivariate analysis of variance (MANOVA):

```
    dataset = "hcp.515"; % LDA on UCLA + ages on 3 bins for HCP
    x = load('../data.515/all_mats.mat');
    y = load('../data.515/all_behav.mat');
    x=x.all_mats;
    y=y.all_behav;
    g = buildGroup(x,dataset,'none'); % mask=false, Bins
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = 10;
    options.phenotype = phenotype('behav',y);
    options.diagnosis = randi(2,175,1);
    options.control=0;
    options.taskID=0;
    options.control = rand(175,1);
    options.phenotypes= y;
    m = manova(g,options);
    m.run();
    m.evaluate();
```

Then in ```manova.m``` under ```run``` function you can see the set of significan edges in ```this.pval(:,1)```.

## Cross Dataset

Cross dataaset evaluation is also included in our pipeline. All you need is to build two groups: ```g1``` and ```g2``` and specify the model for predict (e.g., ```options.model="cpm"```).

```Matlab
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = 2;
    options.model="cpm"; % options: cpm, rcpm
    options.phenotype1 = phenotype('behav',y1);
    options.phenotype2 = phenotype('behav',y2);
    m = cross(g1,g2,options); % first: train, second: test
    m.run();
    m.evaluate();
```
