# Directed Research Course – Numerics.

## About 
This repository contains the final works **I** have carried out in my Directed Research Course (school is left unnamed). The exact material and topics covered will also be left unaddressed for now. Generally speaking, the research revolved around a system of partial differential equations. 

The Matlab scripts here were used in generating the graphs seen in `pde_analysis.pdf`. In particular, the following scripts

- `evolve_diff.m`
- `evolve_diff_G.m`
- `main.m`

`newton_method.m` was a failed but educationally useful test. The rest were used to get numerical solutions to the non-dimensionalized partial differential equations (pde's) we see introduced in `pde_analysis.pdf`. 

There are several problems with the scripts and work listed above. 

1. Lack of documentation.
2. Lack of context.
3. Lack of presentation. 

These scripts weren't necessarily made to be public. That is the primary reason there is a lack of explanations at each step in code. These scripts were used to get results meant to be shown in a group-meeting and in the final report. The physical system in which these equations are trying to model is something that was explained in the meeting and final report. Vaguely speaking, we are measuring the concentration and effect of a substance ("AB") on the surface of a cell. It can't be helped that the `.pdf` looks like a homework problem. 

## Requirements

- GNU Octave, version 4.2.2

## Directions
`main.m` is the only script that needs to be ran. Scripts `evolve_diff.m` and `evolve_diff_G.m` act like functions for `main.m`. The later script is used to solve the 2nd non-dimensionalized pde while the former is used to solve the 1st non-dimensionalized pde. 

The scripts take a considerable amount of time. It’s really dependent on the initial parameters and the density of time-steps and spatial-steps.  

## Final Comments
I don't plan to update the scripts or `.pdf` in any shape or form. It's left **as is**. 

In the future, I might be able to better explain what's going on, but there is no guarantee. 

## Acknowledgments

A big thanks to my course supervisor Mike to publicly display my work.
