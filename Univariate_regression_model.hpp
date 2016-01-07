#ifndef UNIVARIATE_REGRESSION_MODEL_HPP
#define UNIVARIATE_REGRESSION_MODEL_HPP

#include "probability_model.hpp"

//#define M_PI 3.14159265358979323846264338327950288419716939937510
#define LOG_PI log(M_PI)//1.14472988584940017414342735135305871164729481291531



class ur_model : public probability_model{

  public:

  ur_model(double, double, double, Data<double> *);
  ~ur_model();

  virtual  double log_likelihood_interval(changepoint *, changepoint *, changepoint * = NULL);
  virtual  double calculate_mean(changepoint *, changepoint *, changepoint * = NULL);
  void estimate_variance(){m_estimate_variance=true;}
  virtual void use_random_mean(int seed);

  private:

   Data<double> * m_data_ur;
   const double m_alpha; //variance prior shape parameter,can't be zero
   const double m_gamma;//variance prior scale parameter, can't be zero
   const double m_v; //regressor prior variance parameter
   double m_inv_v;
   double m_likelihood_term;
   double * m_ysum;
   double * m_ysum2;
   unsigned long long int m_data_points;
   bool m_estimate_variance;//if true, report E[sigma^2] rather than E[mu]
  bool m_prior_mean;
};


#endif
