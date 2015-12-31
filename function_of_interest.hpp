#ifndef FUNCTION_OF_INTEREST_HPP
#define FUNCTION_OF_INTEREST_HPP

#include "particle.hpp"
#include "changepoint.hpp"
#include "probability_model.hpp"

class Function_of_Interest{

  public:

  Function_of_Interest(int, double, double,double,bool,bool,bool,bool=0,bool=0,int=0,bool=0,double=0);
  ~Function_of_Interest();


  void calculate_function(double, double, Particle<changepoint> **, long long int, double *,double,double,int,bool normalise=1,probability_model * =NULL);
  void normalise_function(double,int,int,int,int iters=0);
  long double * get_g(){return m_exp_last_changepoint;}
  long double * get_variance_g(){return m_variance_exp_last_changepoint;}
  double * get_intensity(){return m_intensity;}
 long  double * get_prob(){return m_prob_function_of_interest;}
  double ** get_prob_sequential(){return m_prob_last_changepoint_sequential;}
  double ** get_g_sequential(){return m_exp_last_changepoint_sequential;}
  double get_average_distance(){return m_average_distance;}
  double get_min_distance(){return m_min_distance;}
  double get_variance_distance(){return m_variance_distance;}
  double calculate_variance_prob(int);
  void reset_prob();
  void set_start(int st){m_start_of_sample = st;}
  void write_mean_to_file(const string output_filename = "intensity.txt");

  private:

   int m_grid;
   double m_start;
   double m_end;
   double m_prior;
   bool m_calculate_g;
   bool m_calculate_prob;
   bool m_calculate_intensity;
   bool m_online;
   bool m_update;
   bool m_instantaneous;   
   double m_delta;
   bool m_fixed;
   double m_grid_points;
   double * m_prior_expectation_function_of_interest;
   double * m_prior_sd_function_of_interest;
   double ** m_exp_last_changepoint_sequential;
   double ** m_prob_last_changepoint_sequential;
   long double * m_exp_last_changepoint;
   long double * m_variance_exp_last_changepoint;
   long double * m_prob_function_of_interest;
   double * m_intensity;
   bool m_cont;
   long double m_average_distance;
   long double m_min_distance;
   long double m_average_distance_squared;
   double m_variance_distance;
   int m_start_of_sample;
   int m_sample_size;
   

};


    

#endif
