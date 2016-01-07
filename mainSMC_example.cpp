#include "SMC_PP_MCMC_nc.hpp"
#include "Data.hpp"
#include "Poisson_process_model.hpp"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include "argument_options_smc.hpp"
#include "probability_model.hpp"
#include "SNCP.hpp"
#include "function_of_interest.hpp"
#include "Univariate_regression_model.hpp"
using namespace std;


int main(int argc, char *argv[])
{
  ArgumentOptionsSMC o = ArgumentOptionsSMC();
  o.parse(argc,argv);

  
  Data<unsigned long long int> * dataobj_int = NULL;
  Data<double> * dataobj = NULL;
  if (o.m_model == "pregression") {
    dataobj_int = new Data<unsigned long long int>(o.m_datafile,false);
  } else {
    dataobj = new Data<double>(o.m_datafile,false);
  }

  if (o.m_model == "ur" || o.m_model == "pregression") {
    o.m_start=0;
    o.m_end=dataobj? dataobj->get_cols():dataobj_int->get_cols();
    cout << o.m_end << endl;
  }
  probability_model * ppptr = NULL;
  
  cout << "seed " << o.m_seed << endl;

  double variance_cp_prior = 0; //if using a prior on the Poisson process parameter for the changepoints
  bool dovariable = 0; //for doing a variable sample size approach
  unsigned int num_proposal_histgoram_bins = 40000/o.m_num_intervals; //number of proposal histogram bins for proposing changepoints in the the sncp model
  unsigned int number_of_data_processes = 1;
  bool only_do_mcmc = false; //when doing SMC repeatedly do MCMC on the intervals [t_0,t_i]
  bool calculate_online_estimate_number_of_cps = true;
  bool sample_from_prior = false;
  bool estimate_var_in_ur = false;

  if(o.m_model == "poisson"){
    ppptr = new pp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj);
    //static_cast<pp_model*>(ppptr)->use_alternative_gamma_prior();
    if (o.m_importance_sampling) {
      ppptr->use_random_mean(o.m_seed);
      /*if (sample_from_prior) {
	ppptr->use_prior_mean();
`	}*/
    }
   
  } else if (o.m_model == "ur") {
    ppptr = new ur_model(o.m_gamma_prior_1,o.m_gamma_prior_2,o.m_v,dataobj);
    if(estimate_var_in_ur)
      static_cast<ur_model*>(ppptr)->estimate_variance();
  } else if (o.m_model == "pregression") {
    ppptr = new pp_model(dataobj_int,NULL,o.m_gamma_prior_1,o.m_gamma_prior_2);
  } else {
    ppptr = new sncp_model(o.m_gamma_prior_1,o.m_gamma_prior_2,dataobj,o.m_seed);
  }


  if (o.m_model == "sncp") {
    o.m_rejection_sampling = false;
    o.m_spacing_prior = false;
    o.m_importance_sampling = false;
    sample_from_prior = false;
  }

  if (sample_from_prior) {
    o.m_rejection_sampling = true;
    o.m_spacing_prior = false;
  }

  if (o.m_rejection_sampling && !sample_from_prior) {
    o.m_spacing_prior = true;
  }


  SMC_PP_MCMC SMCobj(o.m_start, o.m_end, o.m_num_intervals,o.m_particles,o.m_particles,o.m_cp_prior,variance_cp_prior,(probability_model**)&ppptr,number_of_data_processes,dovariable,o.m_calculate_filtering_mean,calculate_online_estimate_number_of_cps,only_do_mcmc, o.m_rejection_sampling, o.m_seed);

  if (o.m_spacing_prior) {
    SMCobj.use_spacing_prior();
  }

  if (o.m_model == "ur") {
    SMCobj.set_discrete_model();
  }

  if(o.m_importance_sampling && o.m_model != "sncp"){
    SMCobj.do_importance_sampling();
  }

  if(o.m_calculate_filtering_mean){
    SMCobj.initialise_function_of_interest(o.m_grid,0,0);
  }

  if (sample_from_prior) {
    SMCobj.sample_from_prior();
  }

  if(!(o.m_disallow_empty_intervals_between_cps || o.m_model == "sncp")){
    SMCobj.set_neighbouring_intervals(1);
  }


  if(o.m_model == "sncp"){
    SMCobj.non_conjugate();
    SMCobj.set_RJ_parameters(o.m_thinning,o.m_burnin,o.m_move_width,"Histogram",(void*)(&num_proposal_histgoram_bins));
  }else{
    SMCobj.set_RJ_parameters(o.m_thinning,o.m_burnin,o.m_move_width);
  }

  SMCobj.set_look_back(1);

  SMCobj.set_ESS_threshold(o.m_ESS_threshold);
  
  if(o.m_print_ESS && !only_do_mcmc){
    SMCobj.store_ESS();
  }

  SMCobj.run_simulation_SMC_PP();


  if(o.m_calculate_filtering_mean){
    SMCobj.print_intensity(0,"intensitySMC.txt");
  }

  if(calculate_online_estimate_number_of_cps)
    SMCobj.print_size_of_sample(0,"kSMC.txt");
  
  if(calculate_online_estimate_number_of_cps)
    SMCobj.print_last_changepoints(0,"taukSMC.txt");
  

  if(o.m_write_cps_to_file){
    SMCobj.print_sample_A(0);
    SMCobj.print_weights();
    SMCobj.print_size_sample_A(0);
    SMCobj.calculate_function_of_interest(o.m_start,o.m_end);
    SMCobj.print_intensity(0,"finalintensitySMC.txt");
  }

  if (o.m_spacing_prior) {
    SMCobj.print_zero_weights(0, "number_zero_weights.txt");
  }

  if (o.m_rejection_sampling) {
    SMCobj.print_rejection_sampling_acceptance_rates(0, "acceptance_rates.txt");
  }

  if(o.m_print_ESS){
    SMCobj.print_ESS(0,"ess.txt");
  }
  if (dataobj)
    delete dataobj;
  //  if (dataobj_int)
  //    delete dataobj_int;
  delete ppptr;  
  return(0);
}
