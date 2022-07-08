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
1. Input values for the following variables depending on the problem (starting from line 18). 
  * `Nnot`: value of $n_0$.
  * `NumMacro`: number of macro-replications.
  * `NumSys`: number of systems considered.
  * `NumConstraint`: number of constraints of the systems
  * `NumThreshold`: maximum number of thresholds of all constraints.
  
2. Choose the value for IZ-parameter $\epsilon_\ell$ and thresholds $q_{\ell m}$, where 
$\ell=1,\ldots,s$ and $m=1,\ldots,d_\ell$. Modify function `configuration()` (starting from line 366) for the following variables. 
  * `epsilon[j]`: values of $\epsilon_j$.
  * `q[m][j]`: values of $q_{j m_j}$.
  
3. Choose the shape of continuation region and compute $\eta_\ell$, where 
$\ell=1,\ldots,s$ as described in the paper. Modify the values of variable `eta[j]` (starting from line 71). Depending on whether the continuation region is chosen as triangular-shaped or straight-line, modify `R[i][j]` (lines 149 and 195). 

4. Decide the correlation between constraints. Modify the following variables from function `read_chol_matrix()`. 
  * Choose the correct input file `cholMatrix_rhoX.txt`, where X is the chosen correlation between each pair of constraints. 
  * Depending on the number of constraints considered, choose the corrent number of inputs for variable `chol_matrix`.

5. Decide the variance configuration between constraints. Modify variable `system_info[i]` from `configuration()` (line 382). We set `system_info[i]=1,2,3` for CV, IV, DV, respectively.  

## 2. Inventory Example

We provide the folder `inventory` that contains some of the required input and code to test the performance of 

The inventory example requires an additional files that contains the mean values for the two performance measures we considered. 
