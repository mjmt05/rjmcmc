#ifndef PROBABILITY_MODEL_HPP
#define PROBABILITY_MODEL_HPP

#include <gsl/gsl_sf_gamma.h>
#include "gsl/gsl_cdf.h"
#include "changepoint.hpp"
#include "Data.hpp"
#include "step_function.hpp"
#include <algorithm>
#include "particle.hpp"
#define LOG_TWO log(2.0)

class probability_model{

  public:

  probability_model(Data<double> * m_data, Step_Function* seasonal_scale = NULL);
  probability_model(vector<string>* data_filenames, double start = 0, double end = DBL_MAX, double season = DBL_MAX, bool make_time_scale = true );
  probability_model(string* data_filename = NULL, string* seasonal_data_filename = NULL, double season = DBL_MAX, bool make_time_scale = true );
  virtual ~probability_model();
  void construct();
  virtual void set_prior_parameters(changepoint *obj1, changepoint* obj2){}
  virtual double log_likelihood_interval(changepoint *, changepoint *, changepoint * = NULL) = 0 ;
  virtual double log_likelihood_interval(double t1, double t2){ return 0;}
  virtual double log_likelihood_interval_with_count(double t1, double t2, unsigned long long int r){return 0;}
  virtual void propose_new_parameters(Particle<changepoint>*, int, unsigned int,changepoint *, changepoint *){};//if third argument 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double calculate_prior_ratio(Particle<changepoint>*,unsigned int){return 0;};//if argument 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double proposal_ratio(Particle<changepoint>*,unsigned int){return 0;};//if 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double calculate_log_predictive_df(double t1, double t2, double t3, bool lower_tail = true ){ return 0; }
  virtual void calculate_sequential_log_predictive_dfs(double start, double end, double increment, bool lower_tail = true, bool two_sided = false, double control_chart_weight = 0.05, string* filename_ptr = NULL, vector<double>* dfs = NULL ){if(dfs) for(unsigned int i=0; i<(end-start)/increment; i++) dfs->push_back(-LOG_TWO);}
  virtual double calculate_log_predictive_df_bounds( double dt, bool lower_tail = true, bool two_sided = false, bool increment_parameters = true ){ m_pvalue_pair = make_pair(-LOG_TWO,-LOG_TWO); return -LOG_TWO; }
  virtual void calculate_posterior_mean_parameters(changepoint *, changepoint *){}
  virtual double calculate_mean(changepoint *, changepoint *, changepoint * = NULL) = 0;
  virtual double draw_mean_from_posterior(changepoint *, changepoint *, changepoint * = NULL){ return m_mean;}
  virtual void set_data_index(changepoint *, unsigned int=0,changepoint * = NULL, changepoint * = NULL);
  virtual double get_mean_function( double t ){ return 1; }
  virtual void collapse_to_seasons_implementation(){}
  virtual void set_parameters_to_current_t(){}
  virtual bool been_active(){ return m_current_data_index>0; }
  virtual bool currently_observable(){ return m_currently_observable; }
  /*SMC functions*/
  virtual void propose_combined_parameters(Particle<changepoint>*,Particle<changepoint>*,changepoint *, double){};/*int tells the index of the particle for the combined region*/
  virtual double non_conjugate_weight_terms(Particle<changepoint>*){return 0;}
  /*return prior parameters*/
  virtual double get_alpha(){return 0;}
  virtual double get_beta(){return 0;}
  virtual void use_random_mean(int){}
  virtual void use_prior_mean(){}
  double get_start() const {return m_start;}
  double get_end() const {return m_end;}
  double get_mean() const {return m_mean;}
  double get_var() const {return m_var;}
  Data<double> * get_data()const{return m_data_cont;}
  void construct_time_scale(vector<string>* data_filenames, double = DBL_MAX );
  Step_Function* get_seasonal_step_function(){ return m_seasonal_scale;}
  void find_data_seasons();
  void collapse_to_seasons();
  double get_two_sided_discrete_p_value( unsigned int );
  double combine_p_values_from_endpoints();
  double combine_p_values_from_endpoints( bool monte_carlo, unsigned int monte_carlo_size = 1);
  void set_parameters_to_current_t( double t ){ m_current_t = t; set_parameters_to_current_t(); }
  pair<double,double> get_p_value_bounds(){ return m_pvalue_pair; }
  bool get_p_value_bounds_on_log_scale(){ return m_pvalue_pair_on_log_scale; }
  static void do_seasonal_analysis(){ m_seasonal_analysis = true; }
  static void do_alternative_p_values(){ m_p_value_alternative_style = true; }
  static void use_monte_carlo_p_values( unsigned int num_samples ){ m_monte_carlo_p_values = true; m_p_value_monte_carlo_samples = num_samples;}
  virtual double log_likelihood_changepoints( unsigned int process, vector<unsigned long long int>& regime_changepoints_data_indices, vector<double>& regime_changepoints_changepoint_positions ) { return 0; }
  virtual double log_likelihood_changepoints( vector<unsigned long long int>& regime_changepoints_data_indices, vector<double>& regime_changepoints_changepoint_positions ) { return 0; }
  void sample_segment_means(Particle<changepoint>*);
  bool m_random_mean;
  void read_in_windows(const std::string& windows_filename="windows.txt",const std::string& window_probs_filename = "window_probs.txt");

 protected:
  double m_mean;
  double m_var;
  double m_start;
  double m_end;
  double m_season;
  unsigned long long int m_current_data_index;
  static double m_current_t;
  Data<double> * m_data_cont;
  Step_Function* m_time_scale;
  Step_Function* m_seasonal_scale;
  unsigned int m_num_data_seasons;
  unsigned int* m_data_seasons;
  double m_log_predictive_pdf;
  double m_log_pdf_const;
  double m_log_predictive_df;
  double m_log_predictive_df2;
  double m_predictive_two_sided_df;
  double m_predictive_two_sided_df2;
  double m_minimum_tail;
  vector<double> m_df_values;
  vector<double> m_pdf_values;
  pair<double,double> m_pvalue_pair;//the upper and lower limits of the current p-value (equal for cts RVs))
  bool m_pvalue_pair_on_log_scale;
  vector< pair<double,double> > m_p_value_endpoints;
  vector<bool> m_p_value_endpoints_log_scale;
  double m_survivor_midpoint;
  bool m_currently_observable;
  static bool m_seasonal_analysis;
  static bool m_p_value_alternative_style;//whether to look at event waiting times, rather than counts, if relevant.
  static bool m_monte_carlo_p_values;//if true, then by default use monte carlo averagigng 
  static unsigned int m_p_value_monte_carlo_samples;
  gsl_rng* m_rng;
  bool m_owner_of_data;
  bool m_owner_of_seasonal_scale;
  bool m_owner_of_time_scale;
  unsigned int m_num_windows;
  unsigned int* m_windows;
  double** m_windowed_lhd_contributions;
  double* m_window_mixture_probs;
  double m_window_mixture_probs_sum;
};





#endif
