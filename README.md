# rjmcmc

DESCRIPTION
___________
rjmcmc is C++ code for both a reversible jump Markov chain Monte Carlo (RJMCMC) algorithm and a sequential Monte Carlo algorithm (SMC) to sample changepoints in continuous time processes. Currently available to specifically run on data modelled as a either a poisson process or a shot noise cox process although easily extended to other processes.

VERSION
_______
Version 1.0 uploaded on 09th February 2015.

SOFTWARE REQUIREMENTS
_____________________
Unix and GSL - Gnu scientific library

INSTALLATION AND BUILDING
_________________________
Download source code and run make mainRJ_example or make mainSMC_example in the source directory for use of the RJMCMC code or SMC code respectively. In the source directory ./mainRJ_example -h or ./mainSMC_example -h will describe the options for input and output for running the code. 


DATA FORMAT
___________
The data file should contain space delimited values, refer to the two example data files shot_noise.txt and coal_data_renormalised.txt.

DOCUMENTATION
_____________
Some detailed instructions available in documentation.pdf but mostly still under development, contact authors with any inquiries. 

REFERENCE
_________
This is the source code used in the following paper:

Turcotte, M. J. M and Heard, N. A. Adaptive sequential Monte Carlo for multiple changepoint analysis. 2015
Turcotte, M.

LICENSING
_________
Copyright (c) 2015 Turcotte, Melissa and Heard, Nicholas

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

