#ifndef POISSON_PROCESS_MODEL_HPP
#define POISSON_PROCESS_MODEL_HPP

#include "probability_model.hpp"
#include <stdlib.h>

class pp_model : public probability_model{

  public:

  pp_model(double, double, Data<double> *, Step_Function* = NULL, Step_Function* = NULL);
  pp_model(vector<string>*, double, double, double = DBL_MAX, double = DBL_MAX, double = DBL_MAX );
  pp_model(vector<string>*);
  pp_model(double alpha,double beta, double rate, Data<double> * data);//Shot noise constructor
  pp_model(string, double=.1, double=.1);//for Poisson regression
  pp_model(Data<unsigned long long int>*, Data<double>* = NULL, double=.1, double=.1,double* =NULL);//for Poisson regression
  ~pp_model();
  void construct();
  void poisson_regression_construct();
  void construct_empirical_prior();

   virtual double log_likelihood_interval(changepoint *, changepoint *);
   virtual double log_likelihood_interval(double t1, double t2);
   virtual void calculate_posterior_mean_parameters(changepoint *, changepoint *);
   virtual double draw_mean_from_posterior(changepoint *, changepoint *);
   virtual double calculate_log_predictive_df(double t1, double t2, double t3, bool lower_tail = true );
   virtual void calculate_sequential_log_predictive_dfs(double start, double end, double increment, bool lower_tail = true, bool two_sided = false, double control_chart_weight = 0.05, string* filename_ptr = NULL, vector<double>* dfs = NULL );
   virtual void set_parameters_to_current_t();
   virtual double calculate_log_predictive_df_bounds( double increment, bool lower_tail = true, bool two_sided = false, bool increment_parameters = true );
   virtual double get_mean_function( double t ){ return m_shot_noise_rate > 0 ? m_pp_time_scale->function(t) : 1; }
   virtual double log_likelihood_changepoints( vector<unsigned long long int>&, vector<double>& );
   virtual double get_alpha(){return m_alpha;}
   virtual double get_beta(){return m_beta;}
   double log_likelihood_up_to(double t);
   virtual double log_likelihood_interval_with_count(double t1, double t2, unsigned long long int r);
   double log_likelihood_length_and_count(double t, unsigned long long int r);
   double log_likelihood_length_and_count(){ return log_likelihood_length_and_count(m_t,m_r); }
   double poisson_regression_log_likelihood_interval(unsigned long long int i1, unsigned long long int i2);
   double calculate_mean(changepoint *, changepoint *);
   double calculate_log_posterior_predictive_pdf( double t, unsigned long long int r );
   double calculate_log_posterior_predictive_df( double t, unsigned long long int r, bool lower_tail = true );//posterior predictive distribution function for a future count of r in t units of time, given the current posterior parameters m_alpha_star, m_beta_star.
   double calculate_waiting_times_log_predictive_df( double increment, bool lower_tail, bool two_sided, bool increment_parameters );
   double calculate_event_count_log_predictive_df( double increment, bool lower_tail, bool two_sided, bool increment_parameters );
  virtual void use_random_mean(int seed);
  virtual void use_prior_mean(){m_posterior_mean = 0;}
  
  
  private:
    unsigned long long int m_r;
    double m_t;
    double m_alpha_star;
    double m_beta_star;
    double m_likelihood_term;
    double m_likelihood_term_zero;//simplified constant for cancellation when #events=0;
    bool m_poisson_regression;//true if the model is Poisson regression rather than Poisson process
    Data<unsigned long long int>* m_cum_counts;//for Poisson regression.
    double m_alpha; //lambda shape parameter,can't be zero
    double m_beta;//lambda scale parameter, can't be zero
    double* m_cum_intensity_multipliers;//fixed mulitpliers for the intensity in Poisson regression.
    Univariate_Function* m_pp_time_scale;
    double m_shot_noise_rate;
    
    bool m_posterior_mean;
};


#endif
