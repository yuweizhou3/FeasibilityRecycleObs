# Online Supplement of "Finding Feasible Systems for Subjective Constraints Using Recycled Observations"

#### by Yuwei Zhou<sup>1</sup>, Sigrun Andradottir<sup>1</sup>, Seong-Hee Kim<sup>1</sup>, and Chuljin Park<sup>2</sup>

#### <sup>1</sup> H. Milton Stewart School of Industrial and Systems Engineering, Georgia Institute of Technology, Atlanta, GA

#### <sup>2</sup> Department of Industrial Engineering, Hanyang University, Seoul, South Korea

## Table of Contents
- [Overview](#overview)
- [Code for Proposed Procedures](#procedure)
- [Code for Inventory Example](#inventory)

## <a name="overview"/> Overview

This repository contains the code and the results for the experiments in Section 6 of the paper "Finding Feasible Systems for Subjective Constraints Using Recycled Observations". 

## <a name="procedure"/> Code for Proposed Procedures

We provide the folder `procedures` that contains the code for Procedure ${\cal RF}$ and the other four alternative procedures ($\text{Recycle}^{\cal B}$, $\text{Restart}^{\rm prob}$, $\text{Restart}^{\rm sum}$, and $\text{Restart}^{\rm max}$). The files are as follows. 
1. `rf.cpp`: this code is used to test the performance of Procedure ${\cal RF}$. 
2. `recycle.cpp`: this code is used to test the performance of Procedure $\text{Recycle}^{\cal B}$. 
3. `restart_max.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm max}$.
4. `restart_prob.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm prob}$.
5. `restart_sum.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm sum}$.
Note that the experiments in Section 6 are only contain by procedures ${\cal RF}$, $\text{Recycle}^{\cal B}$, and $\text{Restart}^{\rm max}$.

As the inputs for all the files are the same, we only discuss how to run `rf.cpp`. The rest four code follows the same instruction. 
1. Decide the number of systems, number of constraints, and the number of thresholds on each constraint. Input those information to `NumSys`, `NumConstraint` (line 18). If the number of thresholds are the same for all constraints, input this number for `NumThreshold` (line 20). Otherwise, input the largest number of thresholds considered among all constraints for `NumThreshold` and adjust the variable `surviveThreshold[j]` (line 121) to incorporate the actual number of thresholds on each constraint.  
2. Decide the correlation between each pair of constraints. Adjust function `read_chol_matrix()` (line 288) to incorporate the chosen correlation and the number of constraints. This function reads a matrix obtained from Cholesky decomposition and is further used to generate observations.
3. Decide the threshold values for each constraint and input them into function `configuration()` (line 367), i.e., adjust variable `q[1][m][j]`.
4. Compute the implementation parameter $\eta$ as discussed in the paper. Input the value of $\eta$ to the variable `eta[i]` (line 81). Depending on how the user wants to set $\eta$, the variable `eta[i]` can be based on systems or based on constraints. 
5. Adjust continuation region variable `R[i][j]` (lines 141 and 190) depending on whether the feasibility check is based on triangular-shaped continuation region or straight-line continuation region. 

Once the above setup is done, the code is ready to be compiled and run. 

## <a name="inventory"/> Inventory Example

We provide the folder `inventory` that contains some of the required input and code to test the performance of Procedures ${\cal RF}$, $\text{Recycle}^{\cal B}$ and $\text{Restart}^{\rm max}$. 

### Data

The inventory example requires an additional files that contains the mean values for the two performance measures we considered. 
