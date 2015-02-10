# rjmcmc

##Description
rjmcmc is C++ code for both a reversible jump Markov chain Monte Carlo (RJMCMC) algorithm and a sequential Monte Carlo algorithm (SMC) to sample changepoints in continuous time processes. Currently available to specifically run on data modelled as a either a poisson process or a shot noise cox process although easily extended to other processes.

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


##Data Format
The data file should contain space delimited values, refer to the two example data files shot_noise.txt and coal_data_renormalised.txt.

##Documentation
Some brief description of the algorithm and the model provided in the Documentation subfolder but mostly still under development, contact authors with any inquiries. 

##Reference
This is the source code used in the following paper:

Turcotte, M. J. M and Heard, N. A. Adaptive sequential Monte Carlo for multiple changepoint analysis. 2015

##Licensing
MIT License  Turcotte, Melissa and Heard, Nicholas
