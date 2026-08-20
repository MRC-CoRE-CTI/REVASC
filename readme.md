# REVASC : code for HIGH RISK REVASC trial design

Ian White, Wenyue Li, Becky Turner
20/08/2026

Simulation to calculate sample size for HIGH RISK REVASC trial design

Main file: 
       simrun.R - calls manyreps via wrapper function, and summarises results as error rates and average powers

Calls: 
       gendata.R - generates data for 7 subgroups
       prepdata.R - fits Cox model to each subgroup and returns logHRs and SEs
       anadata.R - takes subgroup-specific logHRs and SEs and returns posterior probs of theta_k<0, using bayesmeta
       manyreps.R - runs the simulation for 1 DGM over many repetitions

Proposed analysis uses files prepdata.R then anadata.R
