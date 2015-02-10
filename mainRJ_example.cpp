#include "RJMCMC_PP.h"
#include "Data.h"
#include "Poisson_process_model.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "argument_options.h"
#include "probability_model.h"
#include "SNCP.h"
using namespace std;


int main(int argc, char *argv[])
{
  ArgumentOptions o = ArgumentOptions();
  o.parse(argc,argv);

  Data<double> * dataobj = new Data<double>(o.m_datafile,false);


  probability_model * ppptr = NULL;
  if(o.m_model == "poisson"){
    ppptr = new pp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj);
  }else{
    ppptr = new sncp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj,o.m_seed);
  }


  int max_cps = 1e9; //if using a different prior for the changepoints theoretically could use have a maximum number of allowed cps in the model, not implemented
  double variance_cp_prior = 0; //if using a prior on the Poisson process parameter for the changepoints
  bool discrete = 0; //is it a discrete or cts time model
  bool dovariable = 0; //for doing a variable sample size approach
  Particle<changepoint> * initialsample = NULL; //used if you want to start the RJMCMC algorithm with an initial set of changepoints
  bool store_sample = 0; //could store the sample when sampling not recommended.
  unsigned int num_proposal_histgoram_bins = 40000; //number of proposal histogram bins for proposing changepoints in the the sncp model

  rj_pp rjpobject(o.m_start,o.m_end,o.m_iterations,max_cps,o.m_move_width,o.m_cp_prior,variance_cp_prior,ppptr,o.m_thinning,o.m_burnin,discrete,dovariable,initialsample,o.m_seed,store_sample);

  if(o.m_disallow_empty_intervals_between_cps || o.m_model == "sncp"){
    rjpobject.disallow_neighbouring_empty_intervals();
  }

  if(o.m_importance_sampling && o.m_model != "sncp"){
    rjpobject.do_importance_sampling();
  }

  if(o.m_model == "sncp"){
    rjpobject.non_conjugate();
    rjpobject.proposal_type("Histogram",(void*)(&num_proposal_histgoram_bins));
  }

  if(o.m_calculate_posterior_mean){
    rjpobject.calculate_intensity();
    rjpobject.initialise_function_of_interest(o.m_grid);
  }

  if(o.m_write_cps_to_file){
    rjpobject.start_printing_sample();
  }
  
  if(o.m_write_histograms_to_file){
    rjpobject.calculate_sample_histogram(false,o.m_grid,true);
  }

  rjpobject.runsimulation();

  rjpobject.print_acceptance_rates();

  cout << endl;
  cout << "MAP dimension: " << rjpobject.get_MAP_dimension() << endl;
  cout << "MAP changepoints and mean of the MAP dimension: " << *rjpobject.get_MAP_dimension_MAP() << endl;

  if(o.m_calculate_posterior_mean){
    Function_of_Interest * foi = rjpobject.get_foi();
    double normalise = rjpobject.get_sum_thinned_importance_weights();
    foi->normalise_function(normalise,0,o.m_grid,0);
    foi->write_mean_to_file();
  }
  
  if(o.m_write_cps_to_file){
    rjpobject.stop_printing_sample();
  }

  if(o.m_write_histograms_to_file){
    rjpobject.write_frequency_counts_to_file();
    rjpobject.write_1d_histogram_to_file();
  }
 
  delete ppptr;  
  return(0);
}
