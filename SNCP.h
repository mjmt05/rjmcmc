#ifndef SNCP_H
#define SNCP_H

#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "probability_model.h"
#include <stdlib.h>
#include "changepoint.h"
#include "decay_function.h"

class sncp_model : public probability_model{


public:

  sncp_model(double, double, Data<double> *, int=1);
  ~sncp_model();

 virtual double log_likelihood_interval(changepoint *, changepoint *);
 virtual double calculate_mean(changepoint *, changepoint *){return 0;}
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

sncp_model::sncp_model(double alpha, double kappa, Data<double> * data, int seed)
:probability_model(data){
  m_alpha = alpha;
  m_kappa = kappa;
  m_inv_kappa=(double)1/m_kappa;
  m_seed=seed;
  m_prior_ratio_log_term = log(m_alpha);
  gsl_rng_env_setup();
  r_type=gsl_rng_default;
  m_r = gsl_rng_alloc(r_type);
  gsl_rng_set(m_r,seed);
  m_pp_time_scale = new Decay_Function(m_kappa);


}

sncp_model::~sncp_model(){
  gsl_rng_free(m_r);
  delete m_pp_time_scale;

}

double sncp_model::log_likelihood_interval(changepoint *obj1, changepoint *obj2){

  
  unsigned long long int i1 = obj1->getdataindex();
  unsigned long long int i2 = obj2->getdataindex();


  unsigned long long int r = i2-i1;

  double t1 = obj1->getchangepoint();
  double t2 = obj2->getchangepoint();
  double length = t2-t1;

  double theta = obj1->getmeanvalue();

  long double likelihood=1-exp(-m_kappa*length);

  likelihood*=m_inv_kappa;

  likelihood*=-theta;

  likelihood+=r*log(theta);

  likelihood+=m_kappa*r*t1;

  return likelihood;

}

void sncp_model::propose_combined_parameters(Particle<changepoint>* A_particle,Particle<changepoint>* B_particle, changepoint * eoi,  double d){

  int dim_A = A_particle->get_dim_theta();
  int dim_B = B_particle->get_dim_theta();
  changepoint * c1;
  changepoint * c2;


  c1 = A_particle->get_theta_component(dim_A-1);
  c2 = B_particle->get_theta_component(-1);
  double theta1 = c1->getmeanvalue();
  double theta2 = c2->getmeanvalue();

  double discontinuity =theta1*exp(-m_kappa*(d-c1->getchangepoint()))-theta2*exp(-m_kappa*(d-c2->getchangepoint()));  
  m_prior_ratio=-m_prior_ratio_log_term-m_alpha*theta2;

  if(dim_B==0){
    c2=eoi;
  }else{
    c2 = B_particle->get_theta_component(0);
  }

  c1=B_particle->get_theta_component(-1);
  m_proposal_ratio=0;
  m_proposal_ratio+= gamma_distribution_calculations(NULL,c1,c2,NULL,0);
 
  double original_mean;
  double new_mean;
  double exp_term;
  
  if(dim_B>0){
    c1 = B_particle->get_theta_component(0);
    for(int i=1; i<dim_B; i++){
      c2=B_particle->get_theta_component(i);
      m_proposal_ratio-= c1->getlikelihood();
      //cout<<m_proposal_ratio<<endl;
      original_mean = c1->getmeanvalue();
      exp_term=exp(-m_kappa*(c1->getchangepoint()-d));
      new_mean = original_mean+discontinuity*exp_term;
      c1->setmeanvalue(new_mean);
      double likelihood = log_likelihood_interval(c1,c2);
      m_proposal_ratio+=likelihood;
      c1->setlikelihood(likelihood);
      c1=c2;
    }
    original_mean = c1->getmeanvalue();
    m_proposal_ratio-=c1->getlikelihood();
    exp_term=exp(-m_kappa*(c1->getchangepoint()-d));
    new_mean = original_mean+discontinuity*exp_term;
    c1->setmeanvalue(new_mean);
    double likelihood = log_likelihood_interval(c1,eoi);
    m_proposal_ratio+=likelihood;
    c1->setlikelihood(likelihood);
  }
}

double sncp_model::non_conjugate_weight_terms(Particle<changepoint>* p){

  return(m_proposal_ratio+m_prior_ratio);

}

void sncp_model::propose_new_parameters(Particle<changepoint> * parti, int position, unsigned int move_type,changepoint * new_value,changepoint * end_changepoint){
  double mean;

  changepoint * cpobj_right=NULL;
  int p;
  (move_type>0)?p=position+1:p=position;
  (!end_changepoint)?cpobj_right=parti->get_theta_component(p):cpobj_right=end_changepoint;

  changepoint * cpobj_left=NULL;
  (move_type<3)?p=position-1:p=position;
  cpobj_left = parti->get_theta_component(p);

  changepoint * cpobj_double_left=NULL;
  if(p+1>0)
    cpobj_double_left=parti->get_theta_component(p-1);

  changepoint * cpobjact=NULL;
  m_proposal_ratio=0;
  /*prior and proposal before parameters have been changed*/
  if(move_type==0/*birth*/){
    m_prior_ratio=m_prior_ratio_log_term+m_alpha*(cpobj_left->getmeanvalue());
    m_proposal_ratio+=gamma_distribution_calculations(cpobj_double_left,cpobj_left,cpobj_right,NULL,0);
  }else if(move_type==1/*death*/){
    cpobjact = parti->get_theta_component(position);
    mean = cpobj_left->getmeanvalue();
    m_prior_ratio=-m_prior_ratio_log_term;
    m_prior_ratio+=m_alpha*(cpobjact->getmeanvalue()+mean-mean*exp(-m_kappa*(cpobjact->getchangepoint()-cpobj_left->getchangepoint())));
    m_proposal_ratio += gamma_distribution_calculations(cpobj_double_left,cpobj_left,cpobjact,cpobj_right,0);
    m_proposal_ratio+= gamma_distribution_calculations(cpobj_left,cpobjact,cpobj_right,NULL,0);
  }else if(move_type==2/*move changepoint*/){
    mean = cpobj_left->getmeanvalue();
    cpobjact= parti->get_theta_component(position);
    m_prior_ratio=m_alpha*(mean+cpobjact->getmeanvalue()-mean*exp(-m_kappa*(cpobjact->getchangepoint()-cpobj_left->getchangepoint())));
    m_proposal_ratio += gamma_distribution_calculations(cpobj_double_left,cpobj_left,cpobjact,cpobj_right,0);
    m_proposal_ratio += gamma_distribution_calculations(cpobj_left,cpobjact,cpobj_right,NULL,0);
  }else if(move_type==3){
    m_prior_ratio=m_alpha*cpobj_left->getmeanvalue();
    m_proposal_ratio += gamma_distribution_calculations(cpobj_double_left,cpobj_left,cpobj_right,NULL,0);
  }else{cerr<<"no valid move"<<endl; exit(1);}
  


  if(move_type==0 || move_type==2){
    m_proposal_ratio-=gamma_distribution_calculations(cpobj_double_left,cpobj_left,new_value,cpobj_right,1);
    m_proposal_ratio-=gamma_distribution_calculations(cpobj_left,new_value,cpobj_right,NULL,1);
  }else if(move_type==1 || move_type==3)
    m_proposal_ratio-=gamma_distribution_calculations(cpobj_double_left,cpobj_left,cpobj_right,NULL,1);

   


  
  /* prior after parameters have been changed*/
  if(move_type==0){
    mean = cpobj_left->getmeanvalue();
    m_prior_ratio-=m_alpha*(mean+new_value->getmeanvalue()-mean*exp(-m_kappa*(new_value->getchangepoint()-cpobj_left->getchangepoint())));
  }else if(move_type==1){
    m_prior_ratio-=m_alpha*(cpobj_left->getmeanvalue());
  }else if(move_type==2){
    mean = cpobj_left->getmeanvalue();
    m_prior_ratio-=m_alpha*(mean+new_value->getmeanvalue()-mean*exp(-m_kappa*(new_value->getchangepoint()-cpobj_left->getchangepoint())));
  }else if(move_type==3){
    m_prior_ratio-=m_alpha*cpobj_left->getmeanvalue();
  }else{cerr<<"no valid move"<<endl; exit(1);}
 
 

  

}

double sncp_model::gamma_distribution_calculations(changepoint * cpobjleft, changepoint * cpobjact, changepoint * cpobjright, changepoint * cpobjbound, bool propose){


  double t1 = cpobjact->getchangepoint();
  double t3 = cpobjright->getchangepoint();
  double t2=0;
  double length2 = t3-t1;
  double previous_height=0;
  if(cpobjleft){
    t2 = cpobjleft->getchangepoint();
    double length1 = t1-t2;    
    previous_height = cpobjleft->getmeanvalue()*exp(-m_kappa*length1);
  }
  
  long double z = -m_inv_kappa*exp(-m_kappa*length2)+m_inv_kappa+m_alpha;
  unsigned long long int i1 = cpobjact->getdataindex();
  unsigned long long int i2 = cpobjright->getdataindex();
  unsigned long long int r = i2-i1;
  double alpha;
  
  alpha = r+1;
  if(cpobjleft)
    alpha-=previous_height*z;
  if(alpha<1){
    alpha=1; 
  }
 
  double m;
  double bound = 0;
  double length_bound=0;

  if(cpobjbound){
    double t4=cpobjbound->getchangepoint();
    m=cpobjbound->getmeanvalue();
    length_bound = t4-t1;
  }else{
    m=cpobjright->getmeanvalue();
    length_bound=length2;
  }

  if(m>0){
    bound = m*exp(m_kappa*length_bound)-previous_height;
  }
 
 
  long double norm = calculate_normalising_constant(alpha,z,bound);

 
  if(propose){
    double new_value=propose_new_parameters(alpha,z,previous_height,bound, norm);

    if(new_value <0 || previous_height<0){
      cerr<<"SNCP.h: proposing a parameter less than zero"<<endl;
    }
    cpobjact->setmeanvalue(new_value+previous_height);
    return log(calculate_pdf(new_value,alpha,z))-log(norm);
  }else{
    double oldmean = cpobjact->getmeanvalue();
    if(oldmean>0){
      if((oldmean-previous_height)<0){
	cerr<<previous_height<<" "<<t1<<" "<<oldmean<<" "<<t2<<" "<<cpobjleft->getmeanvalue()<<endl;
	exit(1);
      }
      return log(calculate_pdf(oldmean-previous_height,alpha,z))-log(norm);
    }
  }

  return 0;

}


double sncp_model::calculate_pdf(double value, double alpha,double z){
  double pdf =  gsl_ran_gamma_pdf(value,alpha,1.0/(double)z);
  return pdf;
}

double sncp_model::calculate_normalising_constant(double alpha,double z,double bound){
  long double norm =1 ; 
  if(bound>0){
   norm = gsl_cdf_gamma_P(bound,alpha,1.0/(double)z);
  }
  return norm;
}

double sncp_model::propose_new_parameters(double alpha,double z,double  previous_height, double bound,double norm){

  gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
  double new_value=NAN;
  int count =0;
  while(new_value!=new_value){
    double u = gsl_ran_flat(m_r,0,1);
    u*=norm;

    new_value = gsl_cdf_gamma_Pinv(u,alpha,1.0/(double)z); 

    //if(new_value!=new_value)
    //  cout<<u<<" "<<norm<<" "<<bound<<" "<<alpha<<" "<<z<<" "<<previous_height<<endl;
    
    count++;
    if(count>20 and bound>0){
      new_value=(previous_height+bound)/2;
      //cout<<new_value<<endl;
    }
  }
  gsl_set_error_handler(old_handler);
  return new_value;
}
				       


#endif
