#include "Data.hpp"
#include "Poisson_process_model.hpp"
#include "RJMCMC_PP.hpp"
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <map>
#include "SMC_PP_MCMC_nc.hpp"
#include "argument_options_vastdata.hpp"
using namespace std;


int main(int argc, char *argv[])
{

  ArgumentOptionsVast o = ArgumentOptionsVast();
  o.parse(argc,argv);


  /*SMC CAN ONLY TAKE BIAS AND FOI AS DIVERGENCE TYPES*/
 
  Divergence_Type divergence_type = BIAS;
  unsigned int num_runs = 1;
  long long int max_lookup_length = 100000;
  bool SMCMC = 0;
  bool calculate_KL = false;
  bool track_KL = false;

  if (!o.m_sample_sizes.empty()) {
    calculate_KL = true;
  }

  unsigned int divergence_grid = 0;
  if(divergence_type == FOI && !divergence_grid ) {
    divergence_grid = 10;
  }
  
  bool calculate_online_estimate_number_of_cps = true;
  bool calculate_function_of_interest = 0;
  bool calculate_g = 0;
  bool calculate_prob_g = 0;
  bool set_delta=0;
  double delta=2;
  if(calculate_function_of_interest){
    calculate_g = 1;
    calculate_prob_g = 1;
  }

 
 
  //read in data
  string line;
  ifstream myfile(o.m_datafile.c_str());
  istringstream iss(line);  
  vector<string> filenames;

  unsigned int num_of_individuals=0;

  line="";
  myfile >> line;
  while(!myfile.eof()){
    num_of_individuals++;
    filenames.push_back(line);
    line="";
    myfile>>line;
  }

  cerr<<"Number of processes "<<num_of_individuals<<endl;
  if( !o.m_fixed_sample_size ){
    o.m_particles *= num_of_individuals;
  }
  
  cerr<<"Total sample size: "<< o.m_particles << ",   Min. sample size: " << o.m_min_iterations<<endl;
  
  vector<string> f;
  probability_model ** ppptr = new probability_model*[num_of_individuals];
  for(unsigned int i=0; i<num_of_individuals; i++){
    f.erase(f.begin(),f.end());  
    f.push_back(filenames[i]);
    f.push_back("timescale.txt");
    f.push_back("seasonality.txt");
    ppptr[i] = new pp_model(&f,o.m_gamma_prior_1,o.m_gamma_prior_2,o.m_start,o.m_end,1);
  }

  filenames.erase(filenames.begin(),filenames.end()); 

  Histogram_Type<changepoint> combined_histogram(o.m_start,o.m_end,o.m_num_bins,(o.m_end-o.m_start)/o.m_num_bins,true,true,true,0);
  Histogram_Type<changepoint> current_histogram(o.m_start,o.m_end,o.m_num_bins,(o.m_end-o.m_start)/o.m_num_bins,true,true,true,0);
  double sum_entropy=0;
 
  SMC_PP_MCMC * SMCobj = NULL;
  unsigned long long int* sample_sizes=NULL;
  Data<unsigned long long int>* sample_sizes_ptr=NULL;
  if(!o.m_sample_sizes.empty()){
    sample_sizes_ptr = new Data<unsigned long long int>(o.m_sample_sizes);
    sample_sizes = (*sample_sizes_ptr)[0];
    o.m_fixed_sample_size = true;
  }

  for(unsigned int run=1; run<=num_runs; run++){
    unsigned int seed = o.m_seed * run;

    SMCobj = new SMC_PP_MCMC(o.m_start, o.m_end, o.m_num_intervals, o.m_particles, o.m_particles, sample_sizes?&sample_sizes:NULL, o.m_cp_prior, 0, ppptr, num_of_individuals, !o.m_fixed_sample_size, o.m_calculate_filtering_mean, calculate_online_estimate_number_of_cps, SMCMC, 0, seed);

    SMCobj->initialise_function_of_interest(o.m_grid, calculate_g, calculate_prob_g, set_delta, delta, 0);

    if (!o.m_fixed_sample_size) {
      SMCobj->store_sample_sizes();
    }

    SMCobj->set_RJ_parameters(o.m_thinning, o.m_burnin, o.m_move_width);

    if(!o.m_fixed_sample_size){
      SMCobj->set_variable_parameters(divergence_type, o.m_loss_type, o.m_num_bins, o.m_min_iterations, max_lookup_length, divergence_grid);
    }

    if (!o.m_disallow_empty_intervals_between_cps) {
      SMCobj->set_neighbouring_intervals(1);
    }
    SMCobj->set_look_back(1);

    SMCobj->set_ESS_threshold(o.m_ESS_threshold);

    if(o.m_print_ESS && !SMCMC){
      SMCobj->store_ESS();
    }
  
    SMCobj->run_simulation_SMC_PP();

    if(!o.m_fixed_sample_size){
      stringstream out_i_sample_sizes;
      string filename;
      out_i_sample_sizes.str("");
      out_i_sample_sizes << "sample_sizes_";
      if ( num_runs > 1 ) {
	out_i_sample_sizes << run;
      }
      out_i_sample_sizes << ".txt";
      filename=out_i_sample_sizes.str();
      SMCobj->print_variable_sample_sizes(filename.c_str());
    }

    if(calculate_function_of_interest){
      stringstream out_p;
      stringstream out_g;
      out_p << "prob_";
      out_g << "g_";
      if (num_runs > 1) {
	out_p << run;
	out_g << run;
      }
      out_p << ".txt";
      out_g << ".txt";
      SMCobj->print_exp((out_g.str()).c_str());
      SMCobj->print_prob((out_p.str()).c_str());
      out_p.str("");
      out_g.str("");
    }
    
    if (o.m_calculate_filtering_mean) {
      stringstream out_i;
      for(unsigned int ds=0; ds<num_of_individuals; ds++){
	out_i << "intensity_" << ds;
	if (num_runs > 1) {
	  out_i << "_run_" << run;
	}
	out_i << ".txt";
	SMCobj->print_intensity(ds, out_i.str().c_str());
	out_i.str("");
      }
    }

    if (o.m_print_ESS) {
      stringstream out_e;
      for (unsigned int ds = 0; ds < num_of_individuals; ds++) {
	out_e << "ess_" << ds ;
	if (num_runs > 1) {
	  out_e << "_run_" << run;
	}
	out_e << ".txt";
	SMCobj->print_ESS(ds, out_e.str().c_str());
	out_e.str("");
      }
    }

    if(calculate_KL){
      unsigned long long int* num_samples =  SMCobj->get_final_sample_size();
      Particle<changepoint> *** samples = SMCobj->get_sample();
      SMCobj->normalise_weights();
      double ** weights = SMCobj->get_sample_weights();
      for(unsigned long long int i=0; i<num_samples[0]; i++){
	current_histogram.calculate_bin(samples[0][i],samples[0][i]->get_dim_theta());
	current_histogram.increment_bin_counts(NULL,weights[0][i]*num_samples[0]);
      }
      sum_entropy+=current_histogram.get_shannon_entropy();
      combined_histogram.add_histogram(&current_histogram);
      current_histogram.reset();
      if(track_KL)
	cout << combined_histogram.get_shannon_entropy()-sum_entropy/(double)run << endl;//comment out if tracking of JSD not wanted
    }
  

    delete SMCobj;
  }

  for(unsigned int i=0; i<num_of_individuals; i++){
    delete ppptr[i];
  }
  delete [] ppptr;
  
  if(!o.m_fixed_sample_size){
    mc_divergence::delete_lookup_arrays(divergence_type);
  }

  if(sample_sizes_ptr)
    delete sample_sizes_ptr;

  if(calculate_KL && !track_KL){
    double combined_entropy=combined_histogram.get_shannon_entropy();
    double average_entrpopy=sum_entropy/(double)num_runs;
    double jsd = combined_entropy-average_entrpopy;
    cout << jsd << endl;
  }
  
  return(0);
}
