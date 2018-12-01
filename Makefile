CXX=g++ $(INCLUDES)
CXXFLAGS=-Wall -Wno-long-long -pedantic -march=native -O3
INCLUDES=-I/opt/local/include #-I/usr/include/gsl 
LDLIBS=-L/opt/local/lib -lgsl -lgslcblas -lm
OBJS=argument_options.o argument_options_smc.o argument_options_vastdata.o RJMCMC_PP.o SMC_PP_MCMC_nc.o probability_model.o Poisson_process_model.o SNCP.o histogram.o mc_divergence.o changepoint.o function_of_interest.o step_function.o rejection_sampling.o Univariate_regression_model.o 
HEADERS=decay_function.hpp univariate_function.hpp RJMCMC.hpp particle.hpp SMC_PP.hpp Data.hpp histogram_type.hpp 

ifeq ($(DEBUG), 1)
	CXXFLAGS += -DDEBUG -ggdb
endif

.PHONY: clean

all: mainRJ_example mainRJ_seasonal_example mainSMC_example mainSMC_vastdata

mainRJ_example: mainRJ_example.cpp $(OBJS) #$(HEADERS)

mainRJ_seasonal_example: mainRJ_seasonal_example.cpp $(OBJS) #$(HEADERS)

mainSMC_example: mainSMC_example.cpp $(OBJS)#$(HEADERS)

mainSMC_vastdata: mainSMC_vastdata.cpp $(OBJS)

%.o: %.cpp %.hpp

clean:
	rm -f mainRJ_example mainRJ_seasonal_example mainSMC_example mainSMC_vastdata *.o
