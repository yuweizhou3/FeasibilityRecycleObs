# Online Supplement of "Finding Feasible Systems for Subjective Constraints Using Recycled Observations."

#### by Yuwei Zhou<sup>1</sup>, Sigrun Andradottir<sup>1</sup>, Seong-Hee Kim<sup>1</sup>, and Chuljin Park<sup>2</sup>

#### <sup>1</sup> H. Milton Stewart School of Industrial and Systems Engineering, Georgia Institute of Technology, Atlanta, GA

#### <sup>2</sup> Department of Industrial Engineering, Hanyang University, Seoul, South Korea

This repository contains the codes for the experiments in Section 6 of the paper "Finding Feasible Systems for Subjective Constraints Using Recycled Observations".

## 1. Proposed Procedures

We provide the folder `procedures` to that contains the codes for the proposed procedure ${\cal RF}$ and the other four competing procedures. We list the corresponding name for each file as follows.

* `rf.cpp`: this code is used to test the performance of Procedure ${\cal RF}$.
* `recycle.cpp`: this code is used to test the performance of Procedure ${\rm Recycle}^{\cal B}$.
* `restart_max.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm max}$.
* `restart_sum.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm sum}$.
* `restart_prob.cpp`: this code is used to test the performance of Procedure ${\rm Restart}^{\rm prod}$.

We use a 32-bit random number generator, `MRG32k3a()`, to generate independent and identically
uniformly distributed random variates over interval (0,1) which are further used to generate normal random variables. The complete code for this random number generator can be found [here](http://simul.iro.umontreal.ca/rng/MRG32k3a.c).
Note that Section 6 only provides the performance of 
${\rm Restart}^{\rm max}$ which outperforms the other two restart procedures. 

We discuss how to set the inputs to run the codes. Note that the codes for all five procedures use the same structure for user inputs, therefore, we only discuss how to set the inputs for `rf.cpp`.
1. Input values for the following variables depend on the problem (starting from line 18). 
  * `Nnot`: value of $n_0$.
  * `NumMacro`: number of macro-replications.
  * `NumSys`: number of systems considered.
  * `NumConstraint`: number of constraints of the systems
  * `NumThreshold`: maximum number of thresholds of all constraints.
  
2. Choose the value for the IZ-parameter $\epsilon_\ell$ and thresholds $q_{\ell m}$, where 
$\ell=1,\ldots,s$ and $m=1,\ldots,d_\ell$. Modify function `configuration()` (starting from line 366) for the following variables. 
  * `epsilon[j]`: values of $\epsilon_j$.
  * `q[m][j]`: values of $q_{j m_j}$.
  
3. Choose the shape of the continuation region and compute $\eta_\ell$, where 
$\ell=1,\ldots,s$ as described in the paper. Modify the values of variable `eta[j]` (starting from line 71). Depending on the type of a chosen continuation region between triangular-shaped and straight-line ones, modify `R[i][j]` (lines 149 and 195). 

4. Decide the correlation between constraints. Modify the following variables from function `read_chol_matrix()`. 
  * Choose a desired input file `cholMatrix_rhoX.txt`, where X is the correlation between each pair of constraints. 
  * Depending on the number of constraints considered, choose a correct number of inputs for variable `chol_matrix` (comment out unnecessary variables).

5. Decide a variance configuration between constraints. Modify variable `system_info[i]` from function `configuration()` (line 382). We set `system_info[i]=1,2,3` for CV, IV, and DV, respectively. 

The users may choose their own random seeds to generate observations for the feasibility check (by modifying the values of seeds in line 36). In order to retrieve the identical experimental results as in the paper, the users need to use the same random seeds as in the codes.

## 2. Inventory Example

We provide the folder `inventory` that contains the codes for different procedures and a required input file for the $(s,S)$ inventory example in Section 6.4. Note that we include two separate cases depending on whether CRN is applied. 

To run the code, the file `TrueValue.txt`, that contains the estimated mean performance of the two performance measures, is required. The values in this file are obtained through a Markov chain model as discussed in the paper. 

To test the case when CRN is applied, we add a function `generate_demand()` (line 99 of `rf.cpp` as an example) to generate demands that are shared by all systems (which allows the systems to receive demands generated based on the same random seeds). We also need to modify function `generate_one_obs()` (line 343 of `rf.cpp`) for the variable `Demand`.

In order to retrieve the identical experimental results as in the paper, we run the code in parallel by setting different random seeds. More specifically, we run ten simulations in parallel and the number of macro-replications for each simulation is 10,000. The random seeds of each simulation run are set separately as in lines 39-47 of `rf.cpp`. After we receive the results from all ten simulation runs, we take the average for the results in Section 6.4.
