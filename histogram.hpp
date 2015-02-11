#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <map>
#include "Data.h"
#include "mc_divergence.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


using namespace std;

class Histogram
{
  //    friend ostream &operator<< <> ( ostream& , const Histogram<T> &);


public:
  Histogram(double start, double end, unsigned int num_bins, double m_bin_width = -1, bool bounded=true, bool weighted=false, bool one_d=false, unsigned int max_dim=0, string* input_filename=NULL, bool calculate_divergence=false, Divergence_Type divergence_type=BIAS, Loss_Function loss_fn=MINIMAX, unsigned int look_up_length=0,bool estimate_autocorrelation=true);
  Histogram( const Histogram& h );
  ~Histogram();
  void construct(double start, double end, double m_bin_width, unsigned int num_bins, bool bounded, bool weighted, bool one_d, bool calculate_divergence=false, Divergence_Type divergence_type=BIAS, Loss_Function loss_fn=MINIMAX, unsigned int look_up_length=0,bool estimate_autocorrelation=true,unsigned int max_dim=0);
  void reset();//reset histogram counts to zero
  vector<unsigned int>* get_current_bin(){ return &m_current_bin; }
  unsigned long long int get_num_samples(){return !m_weighted? m_samples : (int)m_sum_weights;}
  unsigned long long int increment_bin_counts( vector<unsigned int>* bin = NULL, double weight = 1 );
  void increment_histogram_bin_counts( double weight = 1 );
  void increment_histogram_bin_counts_array( double weight = 1 );
  void increment_1d_bin_counts(double weight = 1);
  void increment_num_bin_repeats(){m_num_bin_repeats++;}
  void add_histogram( Histogram* hist );
  map<vector<unsigned int>,long long unsigned int >* get_bin_counts(){ return &m_histogram_bin_counts; }
  map<vector<unsigned int>,double >* get_bin_weights(){ return &m_histogram_bin_weights; }
  void bin_data( Data<double>* X );
  void write_histogram_to_file(const string output_filename);
  void write_histogram_bin_counts_to_file(const string output_filename);
  void write_histogram_bin_counts_array_to_file(const string output_filename);
  void read_histogram_from_file(const string output_filename, bool increment = false);
  void read_histogram_bin_counts_from_file(const string input_filename);
  void read_histogram_bin_counts_array_from_file(const string input_filename);
  void write_1d_histogram_to_file(const string output_filename);
  void write_marginal_histogram_to_file(const string output_filename, unsigned int dimension=1, unsigned int margin=0);
  void write_true_marginal_histogram_to_file(double(*cdf)(double), const string output_filename);
  void set_divergence_var(double var){ m_mc_divergence->set_var(var); }
  void set_divergence_initial_number_of_bins(){ m_mc_divergence->set_initial_number_of_bins(); }
  static void set_divergence_delta(double delta){ mc_divergence::set_delta(delta); }
  double get_divergence();
  void calculate_bin( double val, bool empty_bin_first=false );
  void push( double val );
  double density( double val );
  double get_shannon_entropy();
  double get_shannon_entropy_from_array();
  double get_true_shannon_entropy( double(*cdf)(double) );
  unsigned long long int get_nonempty_bins(){ return m_nonempty_bins; }
  void normalise_histogram();
  void estimate_repeat_prob();
  double get_autocorrelation_prob(){ return m_autocorrelation_prob; }
  void convert_current_bin_to_int();
  unsigned int sample_bin();
  double sample();
  void calculate_sum_exponentiated_differences();
  unsigned int sample_bin_by_differences();
  double sample_by_differences();
  double sampling_by_differences_log_density( double val );
  void set_seed( unsigned int seed ){ gsl_rng_set(m_r,seed); }
  void track_entropy(){ m_track_entropy = true;}

protected:
  map<vector<unsigned int>,long long unsigned int > m_histogram_bin_counts;
  map<vector<unsigned int>,double > m_histogram_bin_weights;
  vector<unsigned int> m_current_bin;
  map<vector<unsigned int>,long long unsigned int >::iterator m_bin_count_iter;
  map<vector<unsigned int>,double >::iterator m_bin_weight_iter;
  unsigned long long int m_current_bin_int;
  unsigned long long int m_current_base;
  unsigned long long int* m_histogram_bin_counts_array;
  double* m_histogram_bin_weights_array;
  unsigned long long int m_histogram_bin_counts_dim;
  unsigned long long int* m_num_bins_powers;
  double m_start;
  double m_end;
  double m_bin_width;
  bool m_bounded;
  unsigned int m_max_dim;
  unsigned int m_dim_num_bins;
  unsigned int m_dim_num_interior_bins;
  unsigned int m_dim_num_boundary_bins;
  unsigned int m_current_dim_bin;
  unsigned int m_current_dim;//the dimension of the current particle
  unsigned long long int m_nonempty_bins;
  unsigned long long int m_current_bin_count;
  unsigned long long int m_1d_samples;
  unsigned long long int m_samples;
  unsigned long long int* m_1d_histogram_bin_counts;
  double * m_1d_histogram_bin_weights;
  mc_divergence* m_mc_divergence;
  Histogram *m_histogram_odd, *m_histogram_even;
  const Histogram* m_histogram_parent;
  bool m_new_bin;
  double m_independent_repeat_prob, m_autocorrelation_prob;
  bool m_weighted;
  double m_current_bin_weight;
  double m_sum_weights;
  double m_1d_sum_weights;
  bool m_estimate_autocorrelation;
  unsigned long long int m_num_bin_repeats;
  double m_entropy;
  double m_delta_entropy;
  bool m_track_entropy;
  gsl_rng * m_r;
  double m_sum_exponentiated_differences;
  double m_log_sum_exponentiated_differences;
  double m_difference_bin_width;
};

Histogram::Histogram(double start, double end, unsigned int num_bins, double bin_width, bool bounded, bool weighted, bool one_d, unsigned int max_dim, string* input_filename, bool calculate_divergence, Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int look_up_length, bool estimate_autocorrelation){
  m_histogram_parent = NULL;
  if(bin_width<0)
    bin_width = (end-start)/num_bins;
  construct(start,end,bin_width,num_bins,bounded,weighted,one_d,calculate_divergence,divergence_type,loss_fn,look_up_length,estimate_autocorrelation,max_dim);
  if(input_filename)
    read_histogram_from_file(*input_filename);
}

Histogram::Histogram( const Histogram& h ){
  m_histogram_parent = &h;
  if(!h.m_mc_divergence)
    construct(h.m_start,h.m_end,h.m_bin_width,h.m_dim_num_bins,h.m_bounded,h.m_weighted,h.m_1d_histogram_bin_counts!=NULL&&m_1d_histogram_bin_weights!=NULL);
  else
    construct(h.m_start,h.m_end,h.m_bin_width,h.m_dim_num_bins,h.m_bounded,h.m_weighted,h.m_1d_histogram_bin_counts!=NULL&&m_1d_histogram_bin_weights!=NULL,true,h.m_mc_divergence->get_divergence_type(),h.m_mc_divergence->get_loss_function(),0,h.m_estimate_autocorrelation,h.m_max_dim);
}

void Histogram::construct(double start, double end, double bin_width, unsigned int num_bins, bool bounded, bool weighted, bool one_d, bool calculate_divergence, Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int look_up_length, bool estimate_autocorrelation, unsigned int max_dim){
  m_start = start;
  m_end = end;
  m_bin_width = bin_width;
  m_bounded = bounded;
  m_weighted = weighted;
  m_dim_num_interior_bins = num_bins;
  m_dim_num_boundary_bins = m_bounded ? 0:2;
  m_dim_num_bins = m_dim_num_interior_bins + m_dim_num_boundary_bins;
  m_estimate_autocorrelation = estimate_autocorrelation;
  m_max_dim = max_dim;
  m_nonempty_bins = 0;
  m_current_bin_count = 0;
  m_current_bin_weight = 0;
  m_current_dim = 0;
  m_samples = m_1d_samples = 0;
  m_sum_weights = m_1d_sum_weights = 0;
  m_num_bin_repeats = 0;
  m_independent_repeat_prob = 0;
  m_histogram_bin_counts_array = NULL;
  m_histogram_bin_weights_array = NULL;
  m_histogram_bin_counts_dim = 0;
  m_entropy = 0;
  m_delta_entropy = 0;
  m_track_entropy = false;
  m_num_bins_powers = NULL;
  m_sum_exponentiated_differences = 0;
  m_log_sum_exponentiated_differences = 0;
  if(m_max_dim){
    m_num_bins_powers = new unsigned long long int[m_max_dim+1];
    m_num_bins_powers[0]=1;
    for(unsigned int i = 1; i <= m_max_dim; i++)
      m_num_bins_powers[i]=m_num_bins_powers[i-1]*m_dim_num_bins;
    m_histogram_bin_counts_dim = m_num_bins_powers[m_max_dim]+1;
    if(!m_weighted)
      m_histogram_bin_counts_array = new unsigned long long int[m_histogram_bin_counts_dim];
    else
      m_histogram_bin_weights_array = new double[m_histogram_bin_counts_dim];
    for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
      if(!m_weighted)
	m_histogram_bin_counts_array[j] = 0;
      else
	m_histogram_bin_weights_array[j] = 0;
  }
  m_1d_histogram_bin_counts = NULL;
  m_1d_histogram_bin_weights = NULL;
  if( one_d ){
    if(!m_weighted)
      m_1d_histogram_bin_counts = new unsigned long long int[m_dim_num_bins];
    else
      m_1d_histogram_bin_weights = new double[m_dim_num_bins];
    for(unsigned int i = 0; i < m_dim_num_bins; i++)
      if(!m_weighted)
        m_1d_histogram_bin_counts[i] = 0;
      else
        m_1d_histogram_bin_weights[i] = 0;
  }
  if( calculate_divergence )
    m_mc_divergence = new mc_divergence(divergence_type,loss_fn,look_up_length);
  else
    m_mc_divergence = NULL;
  if(!m_histogram_parent && divergence_type == ENTROPY){
    m_histogram_odd = new Histogram( *(Histogram*)this );
    m_histogram_even = new Histogram( *(Histogram*)this );
  }
  else
    m_histogram_odd = m_histogram_even = NULL;
  gsl_rng_env_setup();
  m_r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(m_r,0);
}

void Histogram::reset(){
  m_histogram_bin_counts.clear();
  m_histogram_bin_weights.clear();
  if( m_1d_histogram_bin_counts)
    for(unsigned int i = 0; i < m_dim_num_bins; i++)
      m_1d_histogram_bin_counts[i] = 0;
  if( m_1d_histogram_bin_weights)
    for(unsigned int i = 0; i < m_dim_num_bins; i++)
      m_1d_histogram_bin_weights[i] = 0;
  if(m_mc_divergence)
    m_mc_divergence->reset();
  m_nonempty_bins = 0;
  m_current_bin_count = 0;
  m_samples = m_1d_samples = 0;
  m_sum_weights = m_1d_sum_weights = 0;
  if(m_histogram_odd)
    m_histogram_odd->reset();
  if(m_histogram_even)
    m_histogram_even->reset();
}

Histogram::~Histogram(){
  if(m_1d_histogram_bin_counts)
    delete [] m_1d_histogram_bin_counts;
  if(m_1d_histogram_bin_weights)
    delete [] m_1d_histogram_bin_weights;
  if(m_mc_divergence)
    delete m_mc_divergence;
  if(m_histogram_odd)
    delete m_histogram_odd;
  if(m_histogram_even)
    delete m_histogram_even;
  if(m_histogram_bin_counts_array)
    delete [] m_histogram_bin_counts_array;
  if(m_histogram_bin_weights_array)
    delete [] m_histogram_bin_weights_array;
  if(m_num_bins_powers)
    delete [] m_num_bins_powers;

gsl_rng_free(m_r);
}

unsigned long long int Histogram::increment_bin_counts( vector<unsigned int>* bin, double weight ){
  m_samples++;
  if( bin ){
    m_current_bin = *bin;
    if(m_max_dim)
      convert_current_bin_to_int();
  }
  if(!m_max_dim)
    increment_histogram_bin_counts(weight);
  else
    increment_histogram_bin_counts_array(weight);
  if(m_track_entropy)
    m_entropy += m_delta_entropy;
  if( m_new_bin  )
    m_nonempty_bins++;
  if(m_mc_divergence)
    m_mc_divergence->update_divergence(m_current_bin_count);
  if(m_histogram_odd){
    if(m_samples%2)
      m_histogram_odd->increment_bin_counts(&m_current_bin, weight);
    else
      m_histogram_even->increment_bin_counts(&m_current_bin, weight);
  }
  if(m_estimate_autocorrelation && m_samples>1)
    estimate_repeat_prob();
  return m_current_bin_count;
}

void Histogram::increment_histogram_bin_counts(double weight){
  if(!m_weighted){
    m_bin_count_iter = m_histogram_bin_counts.find(m_current_bin);
    m_new_bin = m_bin_count_iter == m_histogram_bin_counts.end();
  }
  else{
    m_bin_weight_iter = m_histogram_bin_weights.find(m_current_bin);
    m_new_bin = m_bin_weight_iter == m_histogram_bin_weights.end();
    m_sum_weights += weight;
  }
  
  if( !m_new_bin  ){
    if(!m_weighted){
      m_bin_count_iter->second++;
      m_current_bin_count = m_bin_count_iter->second;
      if(m_track_entropy && m_current_bin_count>1)
	m_delta_entropy = -m_current_bin_count*log(m_current_bin_count)+(m_current_bin_count-1)*log(m_current_bin_count-1);
    }
    else{
      if(m_track_entropy){
	m_delta_entropy = -(m_bin_weight_iter->second+weight)*log(m_bin_weight_iter->second+weight);
	if(m_bin_weight_iter->second>0)
	  m_delta_entropy += m_bin_weight_iter->second*log(m_bin_weight_iter->second);
      }
      m_bin_weight_iter->second+= weight;
      m_current_bin_weight = m_bin_weight_iter->second;
    }
  }
  else
  {
    if(!m_weighted){
      m_current_bin_count = m_histogram_bin_counts[m_current_bin] = 1;
      m_delta_entropy = 0;}
    else{
      m_current_bin_weight = m_histogram_bin_weights[m_current_bin] = weight;
      if(m_track_entropy)
	m_delta_entropy = -weight*log(weight);
    }
  }
}

void Histogram::increment_histogram_bin_counts_array(double weight){
  if(!m_weighted){
    if(!m_histogram_bin_counts_array[m_current_bin_int])
      m_new_bin = true;
    m_current_bin_count = ++m_histogram_bin_counts_array[m_current_bin_int];
    if(m_track_entropy)
      m_delta_entropy = m_current_bin_count>1 ? -m_current_bin_count*log(m_current_bin_count)+(m_current_bin_count-1)*log(m_current_bin_count-1) : 0;
  }else{
    if(!m_histogram_bin_weights_array[m_current_bin_int])
      m_new_bin = true;
    if(m_track_entropy){
      m_delta_entropy = -(m_histogram_bin_weights_array[m_current_bin_int] + weight)*log(m_histogram_bin_weights_array[m_current_bin_int] + weight);
      if(m_histogram_bin_weights_array[m_current_bin_int]>0)
	m_delta_entropy = m_histogram_bin_weights_array[m_current_bin_int]*log(m_histogram_bin_weights_array[m_current_bin_int]);
    }
    m_histogram_bin_weights_array[m_current_bin_int] += weight;
    m_current_bin_weight = m_histogram_bin_weights_array[m_current_bin_int];
  }
}

void Histogram::increment_1d_bin_counts(double weight){
  if(!m_1d_histogram_bin_counts && !m_1d_histogram_bin_weights)
    return;
  int previous_bin = -1;
  m_1d_samples++;
  m_1d_sum_weights += weight;
  for(unsigned int k=0; k<m_current_dim; k++){
      if( previous_bin != (int)m_current_bin[k] ){
      if(!m_weighted)
        m_1d_histogram_bin_counts[m_current_bin[k]]++;
      else
        m_1d_histogram_bin_weights[m_current_bin[k]] += weight;
      }
    previous_bin = m_current_bin[k];
  }
}

void Histogram::add_histogram( Histogram* hist ){
  map<vector<unsigned int>,unsigned long long int >::iterator hist_iter;
  for(hist_iter=hist->m_histogram_bin_counts.begin();hist_iter!=hist->m_histogram_bin_counts.end();++hist_iter){
    m_bin_count_iter = m_histogram_bin_counts.find(hist_iter->first);
    if( m_bin_count_iter!= m_histogram_bin_counts.end() )
      m_bin_count_iter->second += hist_iter->second;
    else{
      m_histogram_bin_counts[hist_iter->first] = hist_iter->second;
      m_nonempty_bins++;
    }
  }
  map<vector<unsigned int>,double >::iterator hist_iter2;
  for(hist_iter2=hist->m_histogram_bin_weights.begin();hist_iter2!=hist->m_histogram_bin_weights.end();++hist_iter2){
    m_bin_weight_iter = m_histogram_bin_weights.find(hist_iter2->first);
    if( m_bin_weight_iter!= m_histogram_bin_weights.end() )
      m_bin_weight_iter->second += hist_iter2->second;
    else{
      m_histogram_bin_weights[hist_iter2->first] = hist_iter2->second;
      m_nonempty_bins++;
    }
  }
  m_samples += hist->m_samples;
  m_sum_weights += hist->m_sum_weights;
  if( m_1d_histogram_bin_counts && hist->m_1d_histogram_bin_counts ){
    for(unsigned int i = 0; i < m_dim_num_bins; i++)
      m_1d_histogram_bin_counts[i] += hist->m_1d_histogram_bin_counts[i];
    m_1d_samples += hist->m_1d_samples;
  }
  if( m_1d_histogram_bin_weights && hist->m_1d_histogram_bin_weights ){
    for(unsigned int i = 0; i < m_dim_num_bins; i++)
      m_1d_histogram_bin_weights[i] += hist->m_1d_histogram_bin_weights[i];
    m_1d_sum_weights += hist->m_1d_sum_weights;
  }
  if(m_histogram_bin_counts_array && hist->m_histogram_bin_counts_array)
    for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
      m_histogram_bin_counts_array[j]+=hist->m_histogram_bin_counts_array[j];
  if(m_histogram_bin_weights_array && hist->m_histogram_bin_weights_array)
    for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
      m_histogram_bin_weights_array[j]+=hist->m_histogram_bin_weights_array[j];
  m_track_entropy = false;
}

void Histogram::bin_data( Data<double>* X ){
  unsigned int n = X->get_cols();
  for(unsigned int j = 0; j < n; j++){
    if(!m_bounded ||((*X)[0][j] >= m_start && (*X)[0][j] <= m_end )){
	calculate_bin( (*X)[0][j], true );
	increment_bin_counts();
    }
  }
  calculate_sum_exponentiated_differences();
}

void Histogram::write_histogram_to_file(const string output_filename){
  if(!m_max_dim)
    write_histogram_bin_counts_to_file(output_filename);
  else
    write_histogram_bin_counts_array_to_file(output_filename);
}

void Histogram::write_histogram_bin_counts_to_file(const string output_filename){
  ofstream OutputStream(output_filename.c_str(), ios::out);
  if(!m_weighted)
    for(m_bin_count_iter=m_histogram_bin_counts.begin();m_bin_count_iter!=m_histogram_bin_counts.end();++m_bin_count_iter){
      OutputStream<< m_bin_count_iter->second;
      for(unsigned int i = 0; i < m_bin_count_iter->first.size(); i++)
        OutputStream << " " << (m_bin_count_iter->first)[i];
      OutputStream << endl;
    }
  else
    for(m_bin_weight_iter=m_histogram_bin_weights.begin();m_bin_weight_iter!=m_histogram_bin_weights.end();++m_bin_weight_iter){
      OutputStream<< m_bin_weight_iter->second;
      for(unsigned int i = 0; i < m_bin_weight_iter->first.size(); i++)
        OutputStream << " " << (m_bin_weight_iter->first)[i];
      OutputStream << endl;
    }
}

void Histogram::write_histogram_bin_counts_array_to_file(const string output_filename){
  ofstream OutputStream(output_filename.c_str(), ios::out);
  for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
    if(!m_weighted)
      OutputStream << m_histogram_bin_counts_array[j] << endl;
    else
      OutputStream << m_histogram_bin_weights_array[j] << endl;
}

void Histogram::read_histogram_from_file(const string input_filename, bool increment ){
  if(!increment)
    reset();
  if(!m_max_dim)
    read_histogram_bin_counts_from_file(input_filename);
  else
    read_histogram_bin_counts_array_from_file(input_filename);
}

void Histogram::read_histogram_bin_counts_from_file(const string input_filename){
  ifstream InputStream(input_filename.c_str(), ios::in);
  if(InputStream.is_open()){
    while(InputStream.good()){
      string line;
      getline (InputStream,line);
      if(line.size()>0){
	m_nonempty_bins++;
	istringstream iss(line);
        if(!m_weighted){
          iss >> m_current_bin_count;
          m_samples += m_current_bin_count;}
        else{
          iss >> m_current_bin_weight;
          m_sum_weights += m_current_bin_weight;}
	vector<unsigned int> v;
	unsigned int bin_i;
        int previous_bin = -1;
	while (iss >> bin_i){
	  v.push_back(bin_i);
	  if(m_1d_histogram_bin_counts && previous_bin != (int)bin_i){
	    m_1d_histogram_bin_counts[bin_i] += m_current_bin_count;
            previous_bin = (int)bin_i;
          }
	  if(m_1d_histogram_bin_weights && previous_bin != (int)bin_i){
	    m_1d_histogram_bin_weights[bin_i] += m_current_bin_weight;
            previous_bin = (int)bin_i;
          }
	}
        if(!m_weighted){
          m_histogram_bin_counts[ v ] = m_current_bin_count;
          m_1d_samples += m_current_bin_count;}
        else{
          m_histogram_bin_weights[ v ] = m_current_bin_weight;
          m_1d_sum_weights += m_current_bin_weight;}
      }
    }
    InputStream.close();
  }
}

void Histogram::read_histogram_bin_counts_array_from_file(const string input_filename){
  ifstream InputStream(input_filename.c_str(), ios::in);
  if(InputStream.is_open()){
    for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
    if(!m_weighted)
      InputStream >> m_histogram_bin_counts_array[j];
    else
      InputStream >> m_histogram_bin_weights_array[j];
  }
}

void Histogram::write_1d_histogram_to_file(const string output_filename){
  ofstream OutputStream(output_filename.c_str(), ios::out);
//  OutputStream<<setiosflags(ios::fixed);
//  OutputStream.precision((int)ceil(log10(m_1d_samples)));
  for(unsigned int i = 0; i < m_dim_num_bins; i++ )
    if(!m_weighted)
      OutputStream << (double)m_1d_histogram_bin_counts[i]/m_1d_samples << endl;
  else
      OutputStream << m_1d_histogram_bin_weights[i]/m_1d_sum_weights << endl;
  OutputStream.close();
}

void Histogram::write_marginal_histogram_to_file(const string output_filename, unsigned int dimension, unsigned int margin){
  if(margin>=dimension){
    cerr << "Error: margin of histogram to plot is higher than the dimension" << endl;
    exit(1);
  }
  unsigned int num_appropriate_samples = 0;
  for(m_bin_count_iter=m_histogram_bin_counts.begin();m_bin_count_iter!=m_histogram_bin_counts.end();++m_bin_count_iter){
    if(m_bin_count_iter->first.size()==dimension)
      num_appropriate_samples+=m_bin_count_iter->second;
  }
  if(!num_appropriate_samples)
    return;
  unsigned long long int* counts = new unsigned long long int[ m_dim_num_bins ];
  for(unsigned int i = 0; i < m_dim_num_bins; i++ )
    counts[i]=0;
  for(m_bin_count_iter=m_histogram_bin_counts.begin();m_bin_count_iter!=m_histogram_bin_counts.end();++m_bin_count_iter)
    if(m_bin_count_iter->first.size()==dimension)
      counts[m_bin_count_iter->first[margin]] += m_bin_count_iter->second;
  ofstream OutputStream(output_filename.c_str(), ios::out);
  OutputStream<<setiosflags(ios::fixed);
  OutputStream.precision((int)ceil(log10(num_appropriate_samples)));
  for(unsigned int i = 0; i < m_dim_num_bins; i++ )
    OutputStream << (double)counts[i]/num_appropriate_samples << endl;
  OutputStream.close();
  delete [] counts;
}

void Histogram::write_true_marginal_histogram_to_file(double(*cdf)(double), const string output_filename){
  ofstream OutputStream(output_filename.c_str(), ios::out);
  double F_start = 0;
  if(!m_bounded){
      F_start = cdf(m_start);
      OutputStream << F_start << endl;
  }
  double p1=F_start;
  for(unsigned int i = 1; i<= m_dim_num_interior_bins; i++){
      double p2 = cdf(m_start+i*m_bin_width);
      OutputStream << p2-p1 << endl;
      p1 = p2;
  }
  if(!m_bounded)
      OutputStream << 1-p1 << endl;
  OutputStream.close();
}

double Histogram::get_divergence(){
  double divergence = m_mc_divergence->get_divergence();
  if(m_mc_divergence->get_divergence_type()==ENTROPY)
  {
    double even_odd_divergences = m_histogram_odd->m_mc_divergence->get_divergence()+m_histogram_even->m_mc_divergence->get_divergence();
    divergence -= (m_mc_divergence->get_loss_function()==MINIMAX?.5:.25)*even_odd_divergences;
  }
  return divergence;
}

void Histogram::calculate_bin( double val, bool empty_bin_first ){
  if(empty_bin_first){
    m_current_bin.clear();
    m_current_dim = 0;
  }
  if(m_max_dim && !m_current_dim){
    m_current_base = 1;
    m_current_bin_int = 0;
  }
  if(!m_bounded && val < m_start)
    m_current_dim_bin = 0;
  else if(!m_bounded && val >= m_end)
    m_current_dim_bin = m_dim_num_interior_bins + 1;
  else{
    m_current_dim_bin = static_cast<unsigned int>((val-m_start)/m_bin_width);
    if(!m_bounded)
      m_current_dim_bin++;
  }
  m_current_bin.push_back(m_current_dim_bin);
  m_current_dim++;
  if(m_max_dim){
    if(m_current_dim<=m_max_dim){
      m_current_bin_int += m_current_base * m_current_dim_bin;
      m_current_base *= m_dim_num_bins;
    }else
      m_current_bin_int = m_histogram_bin_counts_dim - 1;
  }
}

void Histogram::push( double val ){
  calculate_bin( val, true );
  increment_bin_counts();
}

double Histogram::density( double val ){
  calculate_bin( val, true );
  return m_histogram_bin_counts[m_current_bin]/(m_samples * m_bin_width);
}

double Histogram::get_shannon_entropy(){
  if(m_track_entropy){
    if(!m_weighted)
      m_entropy = m_entropy/m_samples + log(m_samples);
    else
      m_entropy = m_entropy/m_sum_weights + log(m_sum_weights);
    //    return m_entropy;
  }
  if(m_max_dim>0)
    return get_shannon_entropy_from_array();
  double entropy = 0;
  if(!m_weighted){
    for(m_bin_count_iter=m_histogram_bin_counts.begin();m_bin_count_iter!=m_histogram_bin_counts.end();++m_bin_count_iter)
      if(m_bin_count_iter->second>1)
        entropy -= m_bin_count_iter->second*log(m_bin_count_iter->second);
    entropy = entropy/m_samples + log(m_samples);
  }
  else{
    for(m_bin_weight_iter=m_histogram_bin_weights.begin();m_bin_weight_iter!=m_histogram_bin_weights.end();++m_bin_weight_iter)
      entropy -= m_bin_weight_iter->second*log(m_bin_weight_iter->second);
    entropy = entropy/m_sum_weights + log(m_sum_weights);
  }
  m_entropy = entropy;
  return entropy;
}

double Histogram::get_shannon_entropy_from_array(){
  double entropy = 0;
  for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++){
    double p_j = m_weighted ? m_histogram_bin_weights_array[j] : m_histogram_bin_counts_array[j];
    if(p_j>0)
      entropy -= p_j * log(p_j);
  }
  if(!m_weighted)
    entropy = entropy/m_samples + log(m_samples);
  else
    entropy = entropy/m_sum_weights + log(m_sum_weights);
  return entropy;
}

double Histogram::get_true_shannon_entropy( double(*cdf)(double) ){
  double entropy = 0;
  double F_start = 0;
  if(!m_bounded){
      F_start = cdf(m_start);
      if(F_start>0)
          entropy=-F_start*log(F_start);
  }
  double p1=F_start;
  for(unsigned int i = 1; i<= m_dim_num_interior_bins; i++){
      double p2 = cdf(m_start+i*m_bin_width);
      if(p2>p1)
         entropy -= (p2-p1)*log(p2-p1);
      p1 = p2;
  }
  if(!m_bounded){
      double F_end = p1;
      if(F_end<1)
          entropy-=(1-F_end)*log(1-F_end);
  }
  return entropy;
}

void Histogram::normalise_histogram(){
  if(!m_max_dim){
    if(m_weighted)
      for(m_bin_weight_iter=m_histogram_bin_weights.begin();m_bin_weight_iter!=m_histogram_bin_weights.end();++m_bin_weight_iter)
	m_bin_weight_iter->second/=m_sum_weights;
    else{
      for(m_bin_count_iter=m_histogram_bin_counts.begin();m_bin_count_iter!=m_histogram_bin_counts.end();++m_bin_count_iter)
	m_histogram_bin_weights[m_bin_count_iter->first] = m_bin_count_iter->second/(double)m_samples;
      m_weighted = true;
      m_histogram_bin_counts.clear();
    }
  }else{
    if(!m_weighted)
      m_histogram_bin_weights_array = new double[m_histogram_bin_counts_dim];
    for(unsigned long long int j = 0; j < m_histogram_bin_counts_dim; j++)
      if(m_weighted)
	m_histogram_bin_weights_array[j] /= m_sum_weights;
      else
	m_histogram_bin_weights_array[j] = m_histogram_bin_counts_array[j]/(double)m_samples;
    if(!m_weighted){
      m_weighted = true;
      delete [] m_histogram_bin_counts_array;
      m_histogram_bin_counts_array = NULL;
    }
  }
  m_sum_weights = 1;
}

void Histogram::estimate_repeat_prob(){
    double sample_ratio = (m_samples-1.0)/m_samples;
    m_independent_repeat_prob = m_independent_repeat_prob*sample_ratio*sample_ratio + ((2*m_current_bin_count-1.0)/m_samples)/(double)m_samples;
    m_autocorrelation_prob = (m_num_bin_repeats/(m_samples-1.0) - m_independent_repeat_prob)/(1-m_independent_repeat_prob);
    if(m_autocorrelation_prob<0)
       m_autocorrelation_prob = 0;
}

void Histogram::convert_current_bin_to_int(){
  m_current_bin_int = 0;
  for(unsigned int i = 0; i < m_current_bin.size(); i++)
    m_current_bin_int += m_current_bin[i]*m_num_bins_powers[i];
}

unsigned int Histogram::sample_bin(){
  unsigned int sampled_bin = 0;
  double u = gsl_ran_flat(m_r,0,1) * m_samples;
  double cum_prob = m_histogram_bin_counts[vector<unsigned int>(1,sampled_bin)];
  while( u > cum_prob ){
    sampled_bin++;
    cum_prob += m_histogram_bin_counts[vector<unsigned int>(1,sampled_bin)];
  }
  return sampled_bin;
}

double Histogram::sample(){
  return m_start + (sample_bin()+gsl_ran_flat(m_r,0,1))*m_bin_width;
}

void Histogram::calculate_sum_exponentiated_differences(){
  int left_count = m_histogram_bin_counts[vector<unsigned int>(1,0)];
  int count = m_histogram_bin_counts[vector<unsigned int>(1,1)];
  int max_diff = count - left_count;
  m_log_sum_exponentiated_differences = 1.0;
  for(unsigned int bin = 2; bin < m_dim_num_bins; bin++){
    count = m_histogram_bin_counts[vector<unsigned int>(1,bin)];
    int diff = count-left_count;
    if(diff>max_diff){
      m_log_sum_exponentiated_differences*=exp(diff-max_diff);
      max_diff = diff;
      m_log_sum_exponentiated_differences+=1;
    }
    else
      m_log_sum_exponentiated_differences += exp(diff-max_diff);
    left_count = count;
  }
  m_log_sum_exponentiated_differences = log(m_log_sum_exponentiated_differences) + max_diff;
  m_sum_exponentiated_differences = exp(m_log_sum_exponentiated_differences);
  m_difference_bin_width = (m_bin_width*m_dim_num_bins)/(m_dim_num_bins-1);
}

unsigned int Histogram::sample_bin_by_differences(){
  unsigned int bin = 1;
  int left_count = m_histogram_bin_counts[vector<unsigned int>(1,0)];
  int count = m_histogram_bin_counts[vector<unsigned int>(1,1)];
  double max_diff = count - left_count;
  double log_u = log(gsl_ran_flat(m_r,0,1)) + m_log_sum_exponentiated_differences;
  if( log_u <= max_diff )
    return bin;
  double log_cum_prob = 1.0;
  for(bin = 2; bin < m_dim_num_bins; bin++){
    count = m_histogram_bin_counts[vector<unsigned int>(1,bin)];
    int diff = count-left_count;
    if(diff>max_diff){
      log_cum_prob*=exp(diff-max_diff);
      max_diff = diff;
      log_cum_prob+=1;
    }
    else
      log_cum_prob += exp(diff-max_diff);
    left_count = count;
    if(log_u <= log(log_cum_prob) + max_diff )
      return bin;
  }
  cerr << "Error sampling histogram bin" << endl;
  exit(1);
}

double Histogram::sample_by_differences(){
  return m_start + (sample_bin_by_differences()-1+gsl_ran_flat(m_r,0,1))*m_difference_bin_width;
}

double Histogram::sampling_by_differences_log_density( double val ){
  unsigned int bin = static_cast<unsigned int>((val-m_start)/m_difference_bin_width);
  int count = m_histogram_bin_counts[vector<unsigned int>(1,bin+1)];
  int left_count = m_histogram_bin_counts[vector<unsigned int>(1,bin)];
  return count-left_count-m_log_sum_exponentiated_differences-log(m_difference_bin_width);
}



#endif
