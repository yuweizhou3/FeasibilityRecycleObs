# Online Supplement of "Finding Feasible Systems for Subjective Constraints Using Recycled Observations"

#### by Yuwei Zhou<sup>1</sup>, Sigrun Andradottir<sup>1</sup>, Seong-Hee Kim<sup>1</sup>, and Chuljin Park<sup>2</sup>

#### <sup>1</sup> H. Milton Stewart School of Industrial and Systems Engineering, Georgia Institute of Technology, Atlanta, GA

#### <sup>2</sup> Department of Industrial Engineering, Hanyang University, Seoul, South Korea

This repository contains the code and the results for the experiments in Section 6 of the paper "Finding Feasible Systems for Subjective Constraints Using Recycled Observations".

## 1. Proposed Procedures

We provide the folder `procedures` to contain the code for the proposed procedure ${\cal RF}$ and the other four competing procedures. We list the corresponding name for each file as follows.

* `rf.cpp`: this code is used to test the performance of Procedure ${\cal RF}$.
* `recycle.cpp`: this code is used to test the performance of Procedure ${\rm Recycle}^{\cal B}$.
* `restart_max.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm max}$.
* `restart_sum.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm sum}$.
* `restart_prob.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm prod}$.

Note that Section 6 only provides the performance of ${\rm Restart}^{\rm max}$ among three restart procedures. This is because 
${\rm Restart}^{\rm max}$ donimates the other two restart procedures. 

We discuss how to set the inputs to run the code. Note that the code for all five procedures uses the same structure for user's input, therefore, we only discuss how to set inputs for `rf.cpp`.
1. Input values for the following variables depending on the problem. 
  * `NumSys`
  * `NumConstraint`
  * `Nnot`
  * `NumMacro`
2. Choose the value for IZ-parameter $\epsilon_\ell$, where 
$\ell=1,\ldots,s$.
3. Choose the shape of continuation region and compute $\eta_\ell$, where 
$\ell=1,\ldots,s$ as described in the paper.

## 2. Inventory Example

We provide the folder `inventory` that contains some of the required input and code to test the performance of 

The inventory example requires an additional files that contains the mean values for the two performance measures we considered. 
