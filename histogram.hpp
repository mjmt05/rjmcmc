#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

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
#include "Data.hpp"
#include "mc_divergence.hpp"
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




#endif
