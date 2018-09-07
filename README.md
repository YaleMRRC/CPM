# CPM Project

<center>
<img src="/images/hcp.svg"/>
</center>


This Repository is corrosponding to the Connectome-based Predictive Modeling (CPM). 
The main goal of this algorithm is based around finding the most correlated edges in brain connectome with human behavior measures. 
Based on this, there is a a set of edges positively correlating with measures and a set negatively correlating. 
CPM aims to fit a model which is able to predict human behaviors (e.g., fluid intelligence) for testing subjects. 

# Individualized Parcellation
<center>
<img src="/images/individualized-parcellation.svg"/>
</center>

In classic definition, parcellation is meant to pin some areas of human brain corrosponding to particular functionality. From practical point of view this functionality should translate into higher node activation during a task performance. A traditional view for parcellation is based around a unified group-level network construction. There is good reasons why people use this simple approach; it is comparable among subjects for a prediction measure and also scalable with the number of nodes. But, there is no guarantee to preserve individual differences since a single functional atlas may not perfectly fit for all subjects. In [2] we present a new individualized parcellation which is able to associate a specific point of a parcel (i.e., a voxel) to represent it as a node.

[1] Shen, Xilin, et al. "Using connectome-based predictive modeling to predict individual behavior from brain connectivity." nature protocols 12.3 (2017): 506.

[2] Salehi, Mehraveh, et al. "A submodular approach to create individualized parcellations of the human brain." International Conference on Medical Image Computing and Computer-Assisted Intervention. Springer, Cham, 2017.
