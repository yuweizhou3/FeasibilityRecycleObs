# Online Supplement for Finding Feasible Systems for Subjective Constraints Using Recycled Observations

## Table of Contents
- [Overview](#overview)
- [Code for Proposed Procedures](#procedure)
- [Code for Inventory Example](#inventory)

## <a name="overview"/> Overview

This repository contains the code for the experiments in the paper "Finding Feasible Systems for Subjective Constraints Using Recycled Observations". We provide the instruction on how to run the codes for the proposed procedures and also include a discussion on how we generate data and perform the testing for the inventory example. 

## <a name="procedure"/> Code for Proposed Procedures

We provide the folder `procedures` that contain the code for Procedure ${\cal RF}$ and the other four alternative procedures ($\text{Recycle}^{\cal B}$, $\text{Restart}^{\rm prob}$, $\text{Restart}^{\rm sum}$, and $\text{Restart}^{\rm max}$). The files are as follows. 
1. `rf.cpp`: this code is used to test the performance of Procedure ${\cal RF}$. 
2. `recycle.cpp`: this code is used to test the performance of Procedure $\text{Recycle}^{\cal B}$. 
3. `restart_prob.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm prob}$.
4. `restart_sum.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm sum}$.
5. `restart_max.cpp`: this code is used to test the performance of Procedure $\text{Restart}^{\rm max}$.

As the inputs for all the files are the same, we only discuss how to run `rf.cpp`. The rest four code follows the same instruction. 
1. Decide the number of systems, number of constraints, and the number of thresholds on each constraint. Input those information to `NumSys`, `NumConstraint` (line 18). If the number of thresholds are the same for all constraints, input this number for `NumThreshold` (line 20). Otherwise, input the largest number of thresholds considered among all constraints for `NumThreshold`. 
2. Decide the correlation between each pair of constraints. Adjust function `read_chol_matrix()` (line 288) to incorporate the chosen correlation and the number of constraints. This function reads a matrix obtained from Cholesky decomposition and is further used to generate observations.
