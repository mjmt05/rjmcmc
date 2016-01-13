#ifndef SMC_PP_HPP
#define SMC_PP_HPP

#include "particle.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include "changepoint.hpp"
#include "probability_model.hpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "time.h"
#include "float.h"


template<class T>
class SMC_PP
{

public:

  SMC_PP(double = 0, double =1, unsigned int = 1, long long int = 0, unsigned int = 0, unsigned long long int** = NULL, int=1,bool=0,bool=0,bool=0,int=0);
  virtual ~SMC_PP();

  virtual void sample_particles(double, double) = 0;
  virtual void calculate_weights_join_particles(int,int)=0;
  virtual void resample_particles(double,double,int,const char *,int)=0;
  virtual void ESS_resample_particles(double,int)=0;
  virtual void calculate_function_of_interest(double, double)=0;
  void run_simulation_SMC_PP();
  virtual void delete_samples(int);
  int find_max(double *, unsigned long long int);
  void permute_sample();
  void print_sample_A(int);
  Particle<T>*** get_sample(){ return m_sample_A; }
  unsigned long long int* get_final_sample_size(){ return m_sample_size_A;}
  double** get_sample_weights(){ return m_exp_weights; }
  double * get_sum_weights(){return m_sum_exp_weights;}
  void normalise_weights();
  void print_weights();
  double calculate_ESS(int);
  void set_ESS_threshold(double threshold){m_ESS_percentage=threshold; m_ESS_threshold=m_max_sample_size_A*m_ESS_percentage;}
  void print_size_of_sample(int ds, const char *);
  void print_last_changepoints(int ds, const char *);
  void print_ESS(int ds, const char *);
  void calculate_exp_weights(int);
  void sort( int *, unsigned int);
  void swap( int & , int & );
  unsigned int return_num_ESS(){return m_num_ESS;} 
  void store_ESS();
  void store_sample_sizes();
  void print_sample_birth_times(int);
  void print_size_sample_A(int ds);
  void sample_from_prior() {m_sample_from_prior = true;}
  void do_importance_sampling() {m_importance_sampling = 1;}

protected:

    Particle<T> *** m_sample_A;
    Particle<T> *** m_sample_B;
    Particle<T> *** m_sample_dummy;
    double m_start;
    double m_end;
    double m_change_in_time;
    unsigned int m_num_of_intervals;
    unsigned int m_interval;
    unsigned long long int * m_sample_size_A;
    unsigned long long int m_max_sample_size_A;
    unsigned long long int m_max_sample_size_B;
    unsigned long long int * m_sample_size_B;
    unsigned long long int ** m_sample_sizes;
    bool m_store_sample_sizes;
    double ** m_weights;
    double ** m_exp_weights;
    const gsl_rng_type * r_type;
    gsl_rng * r;
    double m_ESS_threshold;
    double m_ESS_percentage;
    double ** m_ESS;
    bool m_store_ESS;
    double * m_max_weight;
    double * m_sum_exp_weights;
    double * m_sum_squared_exp_weights;
    double ** m_cum_exp_weights;
    int m_num_BF_iterations;
    const char * m_BF_resampling_type;
    double ** m_size_of_sample;
    double ** m_last_changepoint;
    int iters;
    int m_num;
    bool m_variable_B;
    bool m_online_num_changepoints;
    bool m_online_last_changepoint;
    bool MCMC_only;
    //number of data sets
    int * m_process_observed;
    int seed;
    unsigned int m_num_ESS;
  bool  m_sample_from_prior;
  bool m_importance_sampling;
  double log_gamma_pdf(double, double, double);
  probability_model ** m_pm;

};
template<class T>
SMC_PP<T>::SMC_PP(double start, double end, unsigned int num_of_intervals, long long int sizeA, unsigned int sizeB, unsigned long long int** sizes, int num_of_data_sets,bool varyB,bool dochangepoint,bool doMCMC, int s)
  :m_start(start),m_end(end),m_num_of_intervals(num_of_intervals),m_max_sample_size_A(sizeA),m_max_sample_size_B(sizeB),m_sample_sizes(sizes),m_num(num_of_data_sets),m_variable_B(varyB),m_online_num_changepoints(dochangepoint),m_online_last_changepoint(dochangepoint),MCMC_only(doMCMC),seed(s)
{
  m_store_sample_sizes=false;
  m_importance_sampling = 0;
  if(m_sample_sizes){
    m_max_sample_size_A = m_sample_sizes[0][0];
    m_max_sample_size_B = m_sample_sizes[0][1];
  }
  long long int sample_size = m_max_sample_size_A;
  if (varyB) {
    sample_size /= m_num;
  }
  m_sample_from_prior = false;
  m_store_ESS=0;
  if(MCMC_only){
    m_sample_dummy=NULL;
  }
  else{
    m_sample_dummy = new Particle<T>**[m_num];
    for (int i = 0; i < m_num; i++) {
      m_sample_dummy[i] = new Particle<T> *[sample_size];
    }
  }
  m_sample_A = new Particle<T> **[m_num];
  if(!MCMC_only){
    m_sample_B = new Particle<T>**[m_num];
    for(int i=0; i<m_num; i++) {
      m_sample_A[i]= new Particle<T>*[sample_size];
    }
  }else{
    m_sample_B=NULL;
  }
   
  m_exp_weights = new double * [m_num];
  for (int i = 0; i < m_num; i++) {
    m_exp_weights[i] = new double[sample_size];
  }
  m_sum_exp_weights=new double[m_num];
  m_sum_squared_exp_weights=new double[m_num];
  m_num_ESS=0;

  for(unsigned long long int i=0; i<sample_size; i++){
    for(int j=0; j<m_num; j++){
      m_exp_weights[j][i]=1;
    }
  }


  if (MCMC_only){
    m_weights=NULL;
    m_cum_exp_weights=NULL;
  } else {
    m_weights = new double *[m_num];
    m_cum_exp_weights = new double *[m_num];
    
    for (int i = 0; i < m_num; i++) {
      m_weights[i] = new double[sample_size];
      m_cum_exp_weights[i] = new double[sample_size];
    }
    
    for(unsigned long long int i=0; i<sample_size; i++){
      for(int j=0; j<m_num; j++){
	m_weights[j][i] = 0;
      }
    }
  }
 
  m_process_observed = new int[m_num];
  if(!MCMC_only){
    for(int i=0; i<m_num; i++)
      m_process_observed[i]=0;
  }else{
    for(int i=0; i<m_num; i++)
      m_process_observed[i]=1;
  }

  m_max_weight = new double[m_num];
  for(int i=0; i<m_num; i++){
    m_sum_exp_weights[i]=(double)(m_max_sample_size_A);
    m_sum_squared_exp_weights[i]=(double)(m_max_sample_size_A);
    m_max_weight[i]=1;
  }
  
  m_sample_size_A = new unsigned long long int[m_num];
  m_sample_size_B = new unsigned long long int[m_num];

  if(!m_variable_B){
    for( int i=0; i<m_num; i++){
      m_sample_size_B[i]=m_max_sample_size_B;
      m_sample_size_A[i]=m_max_sample_size_A;
    }
  }

  m_change_in_time=(m_end-m_start)/(double)m_num_of_intervals;

  m_ESS_percentage=0.5;
  m_ESS_threshold=m_max_sample_size_A*m_ESS_percentage;
  
  r=NULL;
  gsl_rng_env_setup();
  r_type=gsl_rng_default;
  r = gsl_rng_alloc(r_type);
  gsl_rng_set(r,seed);
  srand(seed);

  m_num_BF_iterations=5;
  m_BF_resampling_type  = "Normal";

  if(m_online_num_changepoints){
    m_size_of_sample = new double * [m_num];
    m_size_of_sample[0] = new double[m_num_of_intervals*m_num];
    for(int i=1; i<m_num; i++)
      m_size_of_sample[i]=m_size_of_sample[i-1]+(m_num_of_intervals);
    for(int i=0; i<m_num; i++)
      for(unsigned int j=0; j<m_num_of_intervals; j++)
	m_size_of_sample[i][j]=0;
  }else{m_size_of_sample=NULL;}
  if(m_online_last_changepoint){
    m_last_changepoint = new double * [m_num];
    m_last_changepoint[0] = new double[m_num_of_intervals*m_num];
    for(int i=1; i<m_num; i++)
      m_last_changepoint[i]=m_last_changepoint[i-1]+(m_num_of_intervals);
    for(int i=0; i<m_num; i++)
      for(unsigned int j=0; j<m_num_of_intervals; j++)
	m_last_changepoint[i][j]=0;
  }else{m_last_changepoint=NULL;}

}
    

template<class T>
SMC_PP<T>::~SMC_PP(){
  long long int sample_size = m_max_sample_size_A;
  long long int i;
  if (m_variable_B) {
    sample_size /= m_num;
  }

  if(m_weights&&!MCMC_only){
    for (i = 0; i < m_num; i++) {
      delete [] m_weights[i];
    }
    delete [] m_weights;
  }

  if(m_cum_exp_weights){
    for (i = 0; i < m_num; i++) {
      delete [] m_cum_exp_weights[i];
    }
    delete []  m_cum_exp_weights;
  }
    
  gsl_rng_free(r);
  for (i = 0; i < m_num; i++) {
    delete [] m_exp_weights[i];
  }
  delete [] m_exp_weights;
  delete [] m_sum_exp_weights;
  delete [] m_sum_squared_exp_weights;
  delete [] m_max_weight;
   
  if (m_sample_dummy&&!MCMC_only){
    delete []  m_sample_dummy;
  }

  if(m_online_num_changepoints){
    delete [] m_size_of_sample[0];
    delete [] m_size_of_sample;
  }
  if(m_online_last_changepoint){
    delete [] m_last_changepoint[0];
    delete [] m_last_changepoint;
  }

  delete [] m_sample_A;
  if(!MCMC_only){
    delete [] m_sample_B;
  }

  delete [] m_sample_size_A;
  delete [] m_sample_size_B;

  delete [] m_process_observed;

  if(m_store_ESS){
    delete [] m_ESS[0];
    delete [] m_ESS;
  }
  if(m_store_sample_sizes){
      delete [] m_sample_sizes[0];
      delete [] m_sample_sizes;
  }
}

template<class T>
void SMC_PP<T>::store_ESS(){
  m_store_ESS=1;
  m_ESS = new double * [m_num];
  m_ESS[0] = new double[m_num*m_num_of_intervals];
  
  for(int i=1; i<m_num; i++)
    m_ESS[i]=m_ESS[i-1] + m_num_of_intervals;

}

template<class T>
void SMC_PP<T>::store_sample_sizes(){
  m_store_sample_sizes=true;
  m_sample_sizes = new unsigned long long int * [m_num];
  m_sample_sizes[0] = new unsigned long long int[m_num*m_num_of_intervals];
  for(int i=1; i<m_num; i++)
    m_sample_sizes[i]=m_sample_sizes[i-1]+m_num_of_intervals;
  for(int i=0; i<m_num; i++){
    for(unsigned int j=0; j<m_num_of_intervals;j++)
      m_sample_sizes[i][j]=0;
  }
}

template<class T>
void SMC_PP<T>::run_simulation_SMC_PP(){

  double * ESS;
  double * BF;
  double * old_sum_weights;
  ESS  = new double[m_num];
  BF = new double[m_num];
  old_sum_weights = new double[m_num];
   
  for (m_interval = 0; m_interval<m_num_of_intervals; m_interval++){
    iters=m_interval;
    for(int ds=0; ds<m_num; ds++){
      old_sum_weights[ds] = log(m_sum_exp_weights[ds]) +m_max_weight[ds];
    }
    if (MCMC_only){           
      sample_particles(m_start,m_start+m_change_in_time*(m_interval+1));
    }
    else{
      sample_particles(m_start+m_change_in_time*m_interval,m_start+m_change_in_time*(m_interval+1));
    }
    if(m_interval>0 && MCMC_only==0 && !m_sample_from_prior){
      permute_sample();
    }
    if (!MCMC_only){
      for(int ds=0; ds<m_num; ds++){
	if(m_process_observed[ds]>0){
	  calculate_weights_join_particles(m_interval,ds);
	  if(m_process_observed[ds]>1){       
	    delete_samples(ds);
	  }
	  for(unsigned int j=0; j<m_sample_size_A[ds]; j++){
	    m_sample_A[ds][j] = m_sample_dummy[ds][j];
	  }
	}
      }
    }
    for(int ds=0; ds<m_num; ds++){
      if (!MCMC_only){
	ESS[ds]=calculate_ESS(ds);
	if(m_store_ESS){
	  m_ESS[ds][m_interval] = ESS[m_interval];
	}
	/*BF[ds] = log(m_sum_exp_weights[ds]) + m_max_weight[ds]-old_sum_weights[ds];
     	  if (BF[ds] < log(0.1) && ESS[ds]>m_ESS_threshold){
	  cout<<"BF "<<ds<<" "<<m_start+m_change_in_time*(i)<<endl;
	  resample_particles(m_start+m_change_in_time*(i),m_start+m_change_in_time*(i+1),m_num_BF_iterations,m_BF_resampling_type,ds);
	  }*/
	m_ESS_threshold=m_sample_size_A[ds]*m_ESS_percentage;

	if (ESS[ds]<m_ESS_threshold){
	  m_num_ESS++;
	  // cout<<"ESS: "<<ds<<" "<<m_interval<<" "<<ESS[ds]<<endl;  
	  ESS_resample_particles(m_start+m_change_in_time*(m_interval+1),ds);
	  ESS[ds]=calculate_ESS(ds);
	  resample_particles(m_start,m_start+m_change_in_time*(m_interval+1),5,"Uniform",ds);
	}
      }
    }

    if (MCMC_only) {
      long long int sample_size = m_max_sample_size_A;
      if (m_variable_B) {
	sample_size /= m_num;
      }
      for(unsigned long long int j=0; j<sample_size; j++){
	for(int ds=0; ds<m_num; ds++){
	  m_exp_weights[ds][j]=1;
	}
      }
    }
    calculate_function_of_interest(m_start+m_change_in_time*(m_interval),m_start+m_change_in_time*(m_interval+1));
    if(m_online_num_changepoints){
      for(int ds=0; ds<m_num; ds++){
	if(m_process_observed[ds]>0){      
	  for(unsigned long long int j=0; j<m_sample_size_A[ds]; j++){
	    m_size_of_sample[ds][m_interval]+=  m_sample_A[ds][j]->get_dim_theta()*m_exp_weights[ds][j];
	  }
	  m_size_of_sample[ds][m_interval]/=m_sum_exp_weights[ds];
	}
      }
    }
    if(m_online_last_changepoint){
      for(int ds=0; ds<m_num; ds++){
	if(m_process_observed[ds]>0){      
	  for(unsigned long long int j=0; j<m_sample_size_A[ds]; j++){
	    m_last_changepoint[ds][m_interval]+=  m_sample_A[ds][j]->get_last_theta_component()->getchangepoint()*m_exp_weights[ds][j];
	  }
	  m_last_changepoint[ds][m_interval]/=m_sum_exp_weights[ds];
	}
      }
    }
  }
  delete [] ESS;
  delete [] BF;
  delete [] old_sum_weights;
}


template<class T>
void SMC_PP<T>::permute_sample(){
  for(int ds=0; ds<m_num; ds++){
    if(m_process_observed[ds]>1)
      gsl_ran_shuffle(r,m_sample_B[ds],m_sample_size_B[ds],sizeof(m_sample_B[ds][0]));
  }
}

template<class T>
double SMC_PP<T>::calculate_ESS(int ds){
    double weights_squared;
    double sum_weights_squared=0;
    double effectivess;

    calculate_exp_weights(ds);

    for (unsigned int long long i=0; i<m_sample_size_A[ds]; i++){   
      weights_squared=m_exp_weights[ds][i]*m_exp_weights[ds][i];
      sum_weights_squared+=weights_squared;
    }
    
    m_sum_squared_exp_weights[ds]=sum_weights_squared;
    effectivess = (m_sum_exp_weights[ds])*(m_sum_exp_weights[ds]/sum_weights_squared);
    return(effectivess);
}



template<class T>
void SMC_PP<T>::calculate_exp_weights(int ds){

  m_sum_exp_weights[ds]=0;
  m_max_weight[ds]=m_weights[ds][find_max(m_weights[ds],m_sample_size_A[ds])];

  for (unsigned long long int i=0; i<m_sample_size_A[ds]; i++){
    if (isinf(m_weights[ds][i])) {
      m_exp_weights[ds][i] = 0;
    } else {
      m_exp_weights[ds][i]=exp(m_weights[ds][i]-m_max_weight[ds]);
    }
    m_sum_exp_weights[ds]+=m_exp_weights[ds][i];
    
    m_cum_exp_weights[ds][i]=m_sum_exp_weights[ds];
  }
}

template<class T>   
void SMC_PP<T>::normalise_weights(){
    for (int i=0; i<m_num; i++){
        for (unsigned long long int j=0; j<m_sample_size_A[i]; j++)
            m_exp_weights[i][j]/=m_sum_exp_weights[i];
        m_sum_exp_weights[i] = 1;
    }
}

template<class T>
void SMC_PP<T>::delete_samples(int ds){

  for(unsigned long long int i=m_sample_size_A[ds]-1; i>0; i--){
    if(m_sample_A[ds][i]!=m_sample_A[ds][i-1]){
      delete m_sample_A[ds][i];
    }
  }
  delete m_sample_A[ds][0];
  
  for(unsigned long long int i=0; i<m_sample_size_B[ds]; i++){
    delete m_sample_B[ds][i];
  }
}

template<class T>
    void SMC_PP<T>::sort(int *arraytosort,unsigned int size)
{

  for (int i=0; i<size-1; i++){
    for (int j=0; j<size-1; j++){
      if( arraytosort[j] > arraytosort[j+1] ){
	swap( arraytosort[j], arraytosort[j+1] );
      }
    }
  }
}

template<class T>
void SMC_PP<T>::swap( int &  element1Ptr, int & element2Ptr )
{
  int hold = element1Ptr;
  element1Ptr = element2Ptr;
  element1Ptr = element2Ptr;
  element2Ptr = hold;
}


template<class T>
int SMC_PP<T>::find_max(double * vec, unsigned long long int size)
{
  int max=0;
  for (int i=1; i<size; i++){
    if (vec[i] > vec[max]){
      max=i;
    }
  }

  return(max);
}


template<class T>
void SMC_PP<T>::print_sample_A(int ds){

    ofstream outfile("sampleSMC.txt",ios::out);

    if(!outfile){
        cerr << "Sample file " << outfile << " could not be opened"<<endl;
	return;
        //exit(1);
    }
    for(unsigned long long int i=0; i<m_sample_size_A[ds]; i++){
      outfile<<*m_sample_A[ds][i];
    }
    outfile.close();
}

template<class T>
void SMC_PP<T>::print_size_sample_A(int ds){

    ofstream outfile("sizesampleSMC.txt",ios::out);

    if(!outfile){
        cerr << "Sample file " << outfile << " could not be opened"<<endl;
	return;
        //exit(1);
    }
    for(unsigned long long int i=0; i<m_sample_size_A[ds]; i++){
      int k = m_sample_A[ds][i]->get_dim_theta();
      outfile << k << ' ';
    }
    outfile << endl;
    outfile.close();
}

template<class T>
void SMC_PP<T>::print_weights(){
 
  ofstream outfile("sample_weights.txt", ios::out);

  if(!outfile){
    cerr<<"Sample weights file " << outfile << " could not be opened"<<endl;
    return;
    //exit(1);
  }

  for(int ds=0; ds<m_num; ds++){
    for(unsigned long long int i=0; i<m_sample_size_A[ds]; i++){
      outfile<<m_exp_weights[ds][i]/m_sum_exp_weights[ds]<<' ';
    }
    outfile<<endl;
  }
  outfile.close();
}


template<class T>
void SMC_PP<T>::print_ESS(int ds, const char * file){
  if(!m_store_ESS){
    cerr << "Can not print ESS to file it has not been stored" << endl;
    return;
  }
  ofstream outfile(file,ios::out);

  if(!outfile){
    cerr<<"ESS file " << outfile << " could not be opened"<<endl;
    return;
    //exit(1);
  }
  for(unsigned int i=0; i<m_num_of_intervals; i++){
    outfile<<m_ESS[ds][i]<<' ';
  }
  outfile<<endl;
  outfile.close();
}

template<class T>
void SMC_PP<T>::print_size_of_sample(int ds, const char * file){
  if(!m_size_of_sample){
    cerr << "Cannot print sample size to file it has not been stored" << endl;
    return;
  }
  
  ofstream outfile(file, ios::out);

  if(!outfile){
    cerr<<"Sample size file " << outfile << " could not be opened"<<endl;
    exit(1);
  }

  for(unsigned int i=0; i<m_num_of_intervals; i++){
    outfile<<m_size_of_sample[ds][i]<<' ';
  }
  outfile<<endl;
  outfile.close();
}

template<class T>
void SMC_PP<T>::print_last_changepoints(int ds, const char * file){
  if(!m_last_changepoint){
    cerr << "Cannot print last changepoint estimates to file, they have not been stored" << endl;
    return;
  }
  ofstream outfile(file, ios::out);
  if(!outfile){
    cerr<<"Last changepoint file " << outfile << " could not be opened"<<endl;
    exit(1);
  }
  for(unsigned int i=0; i<m_num_of_intervals; i++){
    outfile<<m_last_changepoint[ds][i]<<' ';
  }
  outfile<<endl;
  outfile.close();
}



template<class T>
void SMC_PP<T>::print_sample_birth_times(int ds){
  for(unsigned long long int i=m_sample_size_A[ds]-1; i>0; i--){
    if(m_sample_A[ds][i]!=m_sample_A[ds][i-1]){
      cout <<  m_sample_A[ds][i]->get_birth_time() << endl;
    }
  }
 }
    



#endif
