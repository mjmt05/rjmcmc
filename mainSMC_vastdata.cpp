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
//#include "argument_options_vastdata.hpp"
using namespace std;


int main(int argc, char *argv[])
{

  // ArgumentOptionsSMC o = ArgumentOptionsSMC();
  //o.parse(argc,argv);


  /*SMC CAN ONLY TAKE BIAS AND FOI AS DIVERGENCE TYPES*/
 
  Divergence_Type divergence_type = BIAS;
  Loss_Function loss_type = MINIMAX;
  bool calculate_intensity=false;
  bool dovariable= atoi(argv[2]);
  unsigned int num_runs=atoi(argv[3]);
  long long int max_lookup_length = 100000;
  long long int max_iterations = 100000;
  long long int initial_iterations=5000;
  long long int burnin=1000;
  long long int thin=10*5;
  int num_intervals=10;
  int num_bins=275-175-50;//sqrt(max_iterations);//300;
  unsigned int batch_number = 0;
  bool vastdata=1;
  bool store_sample_sizes=1;
  bool SMCMC = 0;
  bool calculate_KL=true;
  bool track_KL=true;
  if(argc>4)
    max_iterations = atol(argv[4]);
  if(argc>5)
    num_bins = atoi(argv[5]);
  if(argc>6)
    divergence_type = static_cast<Divergence_Type>(atoi(argv[6]));
  if(argc>7)
    loss_type = static_cast<Loss_Function>(atoi(argv[7]));
  if(argc>8)
    calculate_intensity=atoi(argv[8]);
  if(argc>9)
    initial_iterations = atoi(argv[9]);
  if(argc>10)
    burnin=atoi(argv[10]);
  if(argc>11)
    thin=atoi(argv[11]);
  if(argc>12)
    batch_number=atoi(argv[12]);
  if(argc>13)
    num_intervals=atoi(argv[13]);
  if(argc>14)
    calculate_KL=atoi(argv[14]);

  
  cerr<<"Number of update interals "<< num_intervals<<endl;
  if(num_runs==1 && dovariable)
    store_sample_sizes=1;


  unsigned int grid_function_of_interest=0;//600;
  if(divergence_type == FOI && !grid_function_of_interest )
    grid_function_of_interest = 100;

  
  if(dovariable && divergence_type == BIAS ){
    cerr<<"Variable sample size bias"<<endl;
  }else if(dovariable && divergence_type == FOI ){
    cerr <<"Variable sample size function of interest (";
    if(calculate_intensity)
      cerr<<"intensity function";
    else
      cerr<<"nearest changepoint function";
    cerr <<")"<<endl;
  }else if(dovariable && divergence_type == CHI_SQUARE ){
    cerr<<"Variable sample size chi squared"<<endl;
  }else if(dovariable && divergence_type == EXTENT ){
    cerr<<"Variable sample size extent"<<endl;
  }else if(dovariable && divergence_type == ENTROPY ){
    cerr<<"Variable sample size JS 50-50"<<endl;
  }else{
    cerr<<"Fixed sample size"<<endl;
  }
  if(loss_type == MINIMAX)
    cerr<<"Minimax Loss" << endl;
  else if(loss_type == AVERAGE)
    cerr<<"Average Loss" << endl;

  //RJ parameters
  double start=0;
  double end=100;
  double nu=0.02;
  double alpha = 1;
  double beta = 0.04;
  double move_width=5;
  end /= 100; move_width /= 500/5; nu *= 100; beta/=8;nu=.1;alpha/=10;beta/=10;
  nu=0.1*10;
  if(vastdata){
    start=0;
    end=10;
    //    nu=1;
    nu=0.005;
    alpha = 0.05;
    beta = 0.05;
    move_width=0.5;
    // num_intervals=240;
  }
  
  bool calculate_function_of_interest=1;
  grid_function_of_interest=num_intervals;
  bool calculate_g=0;
  bool calculate_prob_g=0;

  if(calculate_function_of_interest){
      calculate_g=1;
      calculate_prob_g=1;}

  bool set_delta=0;
  double delta=2;
 
   //read in data
  string line;
  ifstream myfile(argv[1]);
  istringstream iss(line);  
  vector<double> v;
  vector<Data<double>* > vec_data;
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
  //num_of_individuals=5;
  
  if( dovariable ){
    max_iterations *= num_of_individuals;
  }
  //  Histogram combined_histogram_weighted(start,end,num_bins,(end-start)/num_bins,false,false);


  cerr<<"Total sample size: "<<max_iterations<<",   Min. sample size: " <<initial_iterations<<endl;
  cerr<<"Thinning: " << thin << endl;
  cerr<<"Number of runs: "<<num_runs<<endl;
  
cerr<<"Number of bins: "<<num_bins<<endl;


 
  long double * min_sample_sizes = NULL;
  long double * max_sample_sizes = NULL;
  long double * min_sample_sizes_sq = NULL;
  long double * max_sample_sizes_sq = NULL;

  /* if(dovariable){
    min_sample_sizes = new long double[num_of_individuals];
    max_sample_sizes = new long double[num_of_individuals];
    min_sample_sizes_sq = new long double[num_of_individuals];
    max_sample_sizes_sq = new long double[num_of_individuals];

    for(unsigned int i=0; i<num_of_individuals; i++){
      min_sample_sizes[i]=0;  max_sample_sizes[i]=0;  min_sample_sizes_sq[i]=0;  max_sample_sizes_sq[i]=0;
    }
    }*/
 

 
  //double seasonal[3] = {0,0.3314,0.765};
  //double seasonal_mult[3] = {1,5.2,4.1};
  vector<string> f;
  probability_model ** ppptr = new probability_model*[num_of_individuals];
for(unsigned int i=0; i<num_of_individuals; i++){

    f.erase(f.begin(),f.end());  
    f.push_back(filenames[i]);
    f.push_back("timescale.txt");
    f.push_back("seasonality.txt");
    //Data<double> * dataobj = new Data<double>(filenames[i]);
    //ppptr[i] = new pps_model(0.05,0.05,0.005,seasonal,3,seasonal_mult,1,dataobj);
    ppptr[i] = new pp_model(&f,alpha,beta,start,end,1);
 }


if(vastdata)
  filenames.erase(filenames.begin(),filenames.end()); 
 ofstream min_sample_sizes_file;
 ofstream max_sample_sizes_file;
 if(dovariable){
   stringstream out_batch;
   out_batch << "min_sample_sizes_" << batch_number<< ".txt";
   min_sample_sizes_file.open((out_batch.str()).c_str(),ios::out);
   out_batch.str("");
   out_batch << "max_sample_sizes_" << batch_number<< ".txt";
   max_sample_sizes_file.open((out_batch.str()).c_str(),ios::out);
  
 }

 Histogram_Type<changepoint> combined_histogram(start,end,num_bins,(end-start)/num_bins,true,true,true,0);
 Histogram_Type<changepoint> current_histogram(start,end,num_bins,(end-start)/num_bins,true,true,true,0);
 double sum_entropy=0;
 
 SMC_PP_MCMC * SMCobj = NULL;
 unsigned long long int* sample_sizes=NULL;
 Data<unsigned long long int>* sample_sizes_ptr=NULL;
 if(argc>15){
   sample_sizes_ptr = new Data<unsigned long long int>(argv[15]);
   sample_sizes = (*sample_sizes_ptr)[0];
 }

for(unsigned int run=1; run<=num_runs; run++){
  unsigned int runs = run + batch_number*num_runs;
  unsigned int seed = (runs+1) * 1000;

  cerr<<runs<<" "<<num_intervals<<endl;
  SMCobj = new SMC_PP_MCMC(start,end,num_intervals,max_iterations,max_iterations,sample_sizes?&sample_sizes:NULL,nu,0,ppptr,num_of_individuals,dovariable,calculate_intensity,0,SMCMC,0,seed);
  SMCobj->initialise_function_of_interest(grid_function_of_interest,calculate_g,calculate_prob_g,set_delta,delta,0);
  if(store_sample_sizes && dovariable)
    SMCobj->store_sample_sizes();
  SMCobj->set_RJ_parameters(thin,burnin,move_width/(50));
  if(dovariable){
    SMCobj->set_variable_parameters(divergence_type,loss_type,num_bins,initial_iterations,max_lookup_length,grid_function_of_interest);
  }

  SMCobj->set_neighbouring_intervals(1);
  SMCobj->set_look_back(1);
  
  // SMCobj->set_proposal_prior(0.01);
  SMCobj->run_simulation_SMC_PP();

  if(store_sample_sizes && dovariable){
    stringstream out_i_prob;
    string filename;
    out_i_prob.str("");
    out_i_prob<<"sample_sizes_"<<run<<".txt";
    filename=out_i_prob.str();
    SMCobj->print_variable_sample_sizes(filename.c_str());
    
  }

  if(calculate_function_of_interest){

    stringstream out_p;
    stringstream out_g;
    out_p << "prob_" << run << ".txt";
    out_g<<"g_"<<run<<".txt";
    SMCobj->print_exp((out_g.str()).c_str());
    SMCobj->print_prob((out_p.str()).c_str());
    stringstream out_i;
    for(unsigned int ds=0; ds<num_of_individuals; ds++){
      out_i << "intensity_" << ds << ".txt";
      SMCobj->print_intensity(ds, out_i.str().c_str());
      out_i.str("");
    }
    //	SMCobj->print_prob_sequential(ds);
    //SMCobj->calculate_function_of_interest(start,end);
    out_p.str("");
    out_g.str("");
    /*out_p << "prob_final_" << run << ".txt";
      out_g<<"g_final_"<<run<<".txt";
      SMCobj->print_exp((out_g.str()).c_str());
      SMCobj->print_prob((out_p.str()).c_str());*/
  }

  unsigned long long int* num_samples =  SMCobj->get_final_sample_size();
  if(calculate_KL){
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
  
  unsigned long long int * min = NULL;
  if(dovariable) {
    min = SMCobj->get_min_sample_sizes();
    for(unsigned int i=0; i<num_of_individuals; i++){
 
      min_sample_sizes_file<<min[i]<<" ";
      max_sample_sizes_file<<num_samples[i]<<" ";
    }
 
    min_sample_sizes_file<<endl;
    max_sample_sizes_file<<endl;
  }
  delete SMCobj;
 }



cerr<<endl;
 if (dovariable) {
   min_sample_sizes_file.close();
   max_sample_sizes_file.close();
 }



  for(unsigned int i=0; i<num_of_individuals; i++){
    delete ppptr[i];
  }
  delete [] ppptr;
 
  if(dovariable){
 
    delete [] min_sample_sizes;
    delete [] max_sample_sizes;
    delete [] min_sample_sizes_sq;
    delete [] max_sample_sizes_sq;
  }

  if(dovariable){
    mc_divergence::delete_lookup_arrays(divergence_type);
  }
  if(sample_sizes_ptr)
    delete sample_sizes_ptr;

  if(calculate_KL && !track_KL){
    double combined_entropy=combined_histogram.get_shannon_entropy();
    double average_entrpopy=sum_entropy/(double)num_runs;
    double jsd = combined_entropy-average_entrpopy;
    //    cout << combined_entropy << " " << average_entrpopy << " " <<  jsd << endl;
    //    cout << "KL divergence = " << jsd << endl;
    cout << jsd << endl;
  }
  
  return(0);
}
