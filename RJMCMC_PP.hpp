#ifndef RJMCMC_PP_HPP
#define RJMCMC_PP_HPP

#include "RJMCMC.hpp"
#include "changepoint.hpp"
#include "probability_model.hpp"
#include "function_of_interest.hpp"
using namespace std;
#include <string>



using std::ios;


class rj_pp : public rj<changepoint>{

  public:

  rj_pp(double = 0.0 , double =1.0, long long int =10, int = 1000,double = 0, double=1, double=0, probability_model * =NULL, long long int =1,long long int=100,bool=0, bool=0, Particle< changepoint> * = NULL,int=0,bool=false);
  ~rj_pp();

  virtual void initiate_sample(Particle<changepoint> *);
  virtual changepoint* generate_new_parameter()const;
  virtual changepoint* move_parameter(unsigned int) const;
  virtual changepoint * copy_parameter(changepoint *) const;
  virtual double log_likelihood_ratio_birth(changepoint *, int); 
  virtual double log_likelihood_ratio_death(unsigned int);
  virtual double log_likelihood_ratio_move(changepoint *, unsigned int);
  virtual double log_likelihood_ratio_move_parameter(int);
  virtual double log_prior_ratio_birth(changepoint *) const;
  virtual double log_prior_ratio_death(unsigned int) const;
  virtual double log_prior_ratio_move(changepoint *, unsigned int) const;
  virtual double log_prior_ratio_move_parameter() const;
  virtual double log_proposal_ratio_birth(changepoint *) const;
  virtual double log_proposal_ratio_death(int) const;
  virtual double log_proposal_ratio_move() const;
  virtual double log_proposal_ratio_move_parameter() const;
  virtual void update_likelihood_birth(changepoint*,int,double);
  virtual void update_likelihood_death(unsigned int,double);
  virtual void update_likelihood_move(changepoint *, int,double);
  virtual void update_likelihood_move_parameter(int, double);
  virtual void locally_perfect_MAP(Particle<changepoint>* = NULL, bool = false);
  virtual void calculate_function_of_interest();
  void update_primary_function_of_interest();
  void use_one_sided_foi(){m_one_sided_foi = true;}
  Function_of_Interest * get_foi() const {return m_functionofinterest;}
  void disallow_neighbouring_empty_intervals(){m_no_neighbouring_empty_intervals = true;}
  void allow_neighbouring_empty_intervals(){m_no_neighbouring_empty_intervals = false;}
  void calculate_means( Particle< changepoint> * = NULL );
  virtual bool draw_means_from_posterior( bool test_monotonicity = true );
  Step_Function* create_step_function_of_means( Particle< changepoint> * p = NULL, bool normalise = false);
  void write_means_to_file( const char* = "mean_levels.txt", Particle< changepoint> * = NULL, bool normalise = false );//writes the changepoints, the mean levels and the integral of the mean levels to file.
  void write_means_to_ostream( ostream& = std::cout, Particle< changepoint> * = NULL, bool normalise = false );//writes the changepoints, the mean levels and the integral of the mean levels to file.
  void write_changepoints_to_file( ostream&, Particle< changepoint> * = NULL );
  void write_changepoints_and_likelihood_to_file( ostream&, Particle< changepoint> * = NULL );
  virtual void update_new_parameter_with_position( changepoint*, unsigned int );
  virtual double get_component_value( changepoint* cpobj ){ return cpobj->getchangepoint();}
  virtual void set_histogram_component_value_function(){ m_histogram->set_value_function( &changepoint::getchangepoint ); }
  void set_start_cps(double cpstart){m_start_cps=cpstart;log_m_end_time_minus_start_cps=log(m_end_time-m_start_cps);}
  virtual void calculate_importance_weight();
  void calculate_intensity(){m_calculate_mean=1;}
  void initialise_function_of_interest(int,bool=1,bool=1,bool=0,double=0);
  void proposal_type(const char * , void *);
  void set_discrete(bool d){m_discrete=d;}
  void set_spacing_prior() {m_spacing_prior = m_random_nu?false:true;}

  private:
  double m_move_tolerance;
  probability_model * m_pm;
  bool m_calculate_mean;
  double m_prop_distribution_mean;
  double m_prop_distribution_sd;
  changepoint * m_end_of_int_changepoint;
  char m_prop_distribution;
  Histogram * m_prop_histogram;
  void rj_pp_construct();
  Function_of_Interest * m_functionofinterest;
  double log_m_end_time_minus_start_cps;
  bool m_no_neighbouring_empty_intervals;
  bool m_one_sided_foi;
  double m_nu, m_log_nu;//fixed paramater nu for PP prior on changepoints.
  bool m_random_nu;//if true, nu~Gamma(alpha_nu,beta_nu).
  bool m_spacing_prior;
  double m_alpha_nu, m_beta_nu;//gamma hyperparamaters for unknown nu for PP prior on changepoints.
  bool m_discrete;
     
};

#endif
