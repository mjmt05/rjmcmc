#include "Univariate_regression_model.hpp"
#include <gsl/gsl_sf_gamma.h>
#include <math.h>

ur_model::ur_model(double alpha,double gamma, double vconst, Data<double> *data)
:probability_model(),m_alpha(alpha),m_gamma(gamma),m_v(vconst)
{
  m_prior_mean = 0;
    m_data_ur = data;
    m_data_points = m_data_ur->get_cols();
    m_likelihood_term = -0.5*log(m_v)+m_alpha*log(m_gamma)-gsl_sf_lngamma(m_alpha);
    m_inv_v = 1.0 / m_v;
    m_ysum = new double[m_data_points];
    m_ysum2 = new double[m_data_points];
    m_ysum[0]= m_data_ur->m_X[0][0];
    m_ysum2[0] = (m_data_ur->m_X[0][0])*(m_data_ur->m_X[0][0]);
    for (unsigned int i=1; i<m_data_points; i++){
        m_ysum[i]=m_data_ur->m_X[0][i] + m_ysum[i-1];
        m_ysum2[i]=(m_data_ur->m_X[0][i])*(m_data_ur->m_X[0][i]) +m_ysum2[i-1];
    }
    m_estimate_variance=false;
}

ur_model::~ur_model(){
  //    if(m_data_ur)
  //        delete m_data_ur;
    delete [] m_ysum;
    delete[] m_ysum2;
}

void ur_model::use_random_mean(int seed) {
  m_prior_mean = 1;
  m_rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(m_rng, seed);
}

double ur_model::log_likelihood_interval(changepoint *obj1, changepoint *obj2, changepoint *objl1){
    m_current_data_index1=obj1->getdataindex();
    m_current_data_index2=obj2->getdataindex();
    if(m_current_data_index2>m_current_data_index1){
      cerr<<"r can not be less than 0"<<endl;
      exit(1);
    }
    m_r=m_current_data_index2-m_current_data_index1;

    if(m_r==0)
      return 0;
    double like, y, y2;
    like=m_likelihood_term;
    if(m_current_data_index1==0){
      y=m_ysum[m_current_data_index2-1];
      y2=m_ysum2[m_current_data_index2-1];
    }
    else{
      y=m_ysum[m_current_data_index2-1]-m_ysum[m_current_data_index1-1];
      y2=m_ysum2[m_current_data_index2-1]-m_ysum2[m_current_data_index1-1];
    }
    double r_foo = .5*m_r;
    like += -0.5*log(m_inv_v+m_r)-(r_foo)*LOG_PI + gsl_sf_lngamma(r_foo+m_alpha)-(r_foo+m_alpha)*log(m_gamma+0.5*(y2-y*y*(1.0/(m_inv_v+m_r))));
    return(like);
}


double ur_model::calculate_mean(changepoint *obj1, changepoint *obj2, changepoint *objl1){

  if (m_prior_mean) {
    double m_var = gsl_ran_gamma(m_rng, m_alpha, m_gamma);
    m_var = 1.0/m_var;
    m_mean = gsl_ran_gaussian(m_rng, m_var * m_v);
    if (!m_estimate_variance)
      return m_mean;
    return m_var;
  }
  double y, y2, r;
  unsigned long long int dataindex1, dataindex2;
  dataindex1=obj1->getdataindex();
  dataindex2=obj2->getdataindex();
  r = dataindex2-dataindex1;
  if (r<0){
    cerr<<"r can not be less than 0"<<endl;
    exit(1);
  }

  if (r==0){
    m_mean=0;
    m_var=m_gamma/(m_alpha-1);
  }
  else{
    if(dataindex1==0){
      y=m_ysum[dataindex2-1];
      y2=m_ysum2[dataindex2-1];
    }
    else{
      y=m_ysum[dataindex2-1]-m_ysum[dataindex1-1];
      y2=m_ysum2[dataindex2-1]-m_ysum2[dataindex1-1];
    }
    m_mean = (1.0/(m_inv_v+r))*(y);
    m_var = (m_gamma+0.5*(y2-y*y*(1.0/(m_inv_v+r))))/(m_alpha+.5*r-1);
  }
  if(!m_estimate_variance)
    return(m_mean);
  return m_var;
}
