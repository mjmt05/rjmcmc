#include "SMC_PP_MCMC_nc.h"
#include "Data.h"
#include "Poisson_process_model.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "argument_options_smc.h"
#include "probability_model.h"
#include "SNCP.h"
#include "function_of_interest.h"
using namespace std;


int main(int argc, char *argv[])
{
  ArgumentOptionsSMC o = ArgumentOptionsSMC();
  o.parse(argc,argv);

  Data<double> * dataobj = new Data<double>(o.m_datafile,false);


  probability_model * ppptr = NULL;
  if(o.m_model == "poisson"){
    ppptr = new pp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj);
  }else{
    ppptr = new sncp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj,o.m_seed);
  }


  double variance_cp_prior = 0; //if using a prior on the Poisson process parameter for the changepoints
  bool dovariable = 0; //for doing a variable sample size approach
  unsigned int num_proposal_histgoram_bins = 40000/o.m_num_intervals; //number of proposal histogram bins for proposing changepoints in the the sncp model
  unsigned int number_of_data_processes = 1;
  bool only_do_mcmc = 0; //when doing SMC repeatedly do MCMC on the intervals [t_0,t_i]
  bool calculate_online_estimate_number_of_cps = 0;

  SMC_PP_MCMC SMCobj(o.m_start, o.m_end, o.m_num_intervals,o.m_particles,o.m_particles,o.m_cp_prior,variance_cp_prior,(probability_model**)&ppptr,number_of_data_processes,dovariable,o.m_calculate_filtering_mean,calculate_online_estimate_number_of_cps,only_do_mcmc,o.m_seed);

  if(o.m_calculate_filtering_mean){
    SMCobj.initialise_function_of_interest(o.m_grid,0,0);
  }

  if(!(o.m_disallow_empty_intervals_between_cps || o.m_model == "sncp")){
    SMCobj.set_neighbouring_intervals(1);
  }

  if(o.m_importance_sampling && o.m_model != "sncp"){
    Function_of_Interest * foi = SMCobj.get_function_of_interest(0);
    foi->set_importance_sampling();
  }

  if(o.m_model == "sncp"){
    SMCobj.non_conjugate();
    SMCobj.set_RJ_parameters(o.m_thinning,o.m_burnin,o.m_move_width,"Histogram",(void*)(&num_proposal_histgoram_bins));
  }else{
    SMCobj.set_RJ_parameters(o.m_thinning,o.m_burnin,o.m_move_width);
  }
  SMCobj.set_look_back(1);

  SMCobj.set_ESS_threshold(o.m_ESS_threshold);
  
  if(o.m_print_ESS){
    SMCobj.store_ESS();
  }

  SMCobj.run_simulation_SMC_PP();


  if(o.m_calculate_filtering_mean){
    SMCobj.print_intensity(0,"intensitySMC.txt");
  }

  if(o.m_write_cps_to_file){
    SMCobj.print_sample_A(0);
    SMCobj.print_weights();
    SMCobj.print_size_sample_A(0);
    SMCobj.calculate_function_of_interest(o.m_start,o.m_end);
    SMCobj.print_intensity(0,"finalintensitySMC.txt");
  }

  if(o.m_print_ESS){
    SMCobj.print_ESS(0,"ess.txt");
  }
  
 
  delete ppptr;  
  return(0);
}
