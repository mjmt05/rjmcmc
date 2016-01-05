#ifndef SNCP_HPP
#define SNCP_HPP



#include "probability_model.hpp"
#include <stdlib.h>
#include "changepoint.hpp"
#include "particle.hpp"
#include "decay_function.hpp"

class sncp_model : public probability_model{


public:

  sncp_model(double, double, Data<double> *, int=1);
  ~sncp_model();

 virtual double log_likelihood_interval(changepoint *, changepoint *, changepoint * = NULL);
 virtual double calculate_mean(changepoint *, changepoint *, changepoint * = NULL){return 0;}
 virtual void propose_new_parameters(Particle<changepoint>*, int, unsigned int, changepoint *, changepoint *);
 virtual double calculate_prior_ratio(Particle<changepoint>*,unsigned int){return(m_prior_ratio);}
 double gamma_distribution_calculations(changepoint *, changepoint *, changepoint *,changepoint * =NULL, bool=1);
 virtual double proposal_ratio(Particle<changepoint>*,unsigned int){return(m_proposal_ratio);}
 virtual double get_mean_function( double t ){ return m_pp_time_scale->function(t); }
 virtual void propose_combined_parameters(Particle<changepoint>*,Particle<changepoint>*,changepoint *, double);/*int tells the index of the particle for the combined region*/
 virtual double non_conjugate_weight_terms(Particle<changepoint>*);
 double calculate_pdf(double, double, double);
 double calculate_normalising_constant(double , double , double);
 double propose_new_parameters(double, double, double, double, double);
  
private:

 double m_alpha; //prior parameter for the intensity jumps
 double m_kappa; //rate of decline for the intensity between changepoints
 double m_inv_kappa;
 int m_seed;
 double m_prior_ratio;
 double m_prior_ratio_log_term;
 double m_proposal_ratio;
 Decay_Function* m_pp_time_scale;

 const gsl_rng_type * r_type;
 gsl_rng * m_r;
};


#endif
