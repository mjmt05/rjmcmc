# rjmcmc

##Description
rjmcmc is C++ code for both a reversible jump Markov chain Monte Carlo (RJMCMC) algorithm and a sequential Monte Carlo algorithm (SMC) to sample changepoints in continuous time processes. 

##Version
Version 1.0 uploaded on 09th February 2015.

##Software Requirements
Unix and GSL - Gnu scientific library

##Usage
RJMCMC
```
make mainRJ_example
./mainRJ_example -h
```

SMC
```
make mainSMC_example
./mainSMC_example -h
```

VASTDATA 
```
make mainSMC_vastdata
./mainSMC_vastdata -h
```
Use the default examples for SMC to get results from the paper, see reference. To run the vast data copy the executable (mainSMC_vastdata) into the folder vastdata and use the default example for results from the paper.

##Data Format
The data file should contain space delimited values, refer to the two example data files shot_noise.txt and coal_data_renormalised.txt.

##Documentation
Some brief description of the algorithm and the model provided in the Documentation subfolder but mostly still under development, contact authors with any inquiries. The SMC algorithm and most models used are detailed in the paper below.  

##Reference
This is the source code used in the following paper:

Adaptive sequential Monte Carlo for multiple changepoint analysis. 2015 
<!-- Turcotte, M. J. M and Heard, N. A. Adaptive sequential Monte Carlo for multiple changepoint analysis. 2015 -->

Licensing
Copyright &copy;  2015. Melissa Turcotte and Nicholas Heard. Released under the MIT License. 
