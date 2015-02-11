CXX=g++
CXXFLAGS=-Wall -Wno-long-long -pedantic -march=native -O3
INCLUDES=-I/opt/local/include #-I/usr/include/gsl 
LDLIBS=-L/opt/local/lib -lgsl -lgslcblas -lm
OBJS=argument_options.o argument_options_smc.o RJMCMC_PP.o SMC_PP_MCMC_nc.o probability_model.o Poisson_process_model.o SNCP.o histogram.o mc_divergence.o changepoint.o function_of_interest.o step_function.o
HEADERS=decay_function.hpp univariate_function.hpp RJMCMC.hpp particle.hpp SMC_PP.hpp Data.hpp histogram_type.hpp

.PHONY: clean

all: mainRJ_example mainSMC_example

mainRJ_example: mainRJ_example.cpp $(OBJS) $(HEADERS)

mainSMC_example: mainSMC_example.cpp $(OBJS) $(HEADERS)

%.o: %.cpp %.hpp

clean:
	rm -f mainRJ_example mainSMC_example *.o
