#include "RJMCMC_PP.hpp"
#include "RJMCMC.hpp"
#include "string.h"
#include "changepoint.hpp"

rj_pp::rj_pp(double begin, double end, long long int its, int max, double tol, double e_nu, double v_nu, probability_model *ppptr, long long int thin, long long int burnin, bool discrete,bool calc_KL, Particle<changepoint> * initialparticle,int se, bool store_sample)
:rj<changepoint>(begin,end,its,max,thin,burnin,calc_KL,initialparticle,se,store_sample),m_move_tolerance(tol),m_pm(ppptr),m_discrete(discrete)
{   
  m_nu = e_nu;
  m_log_nu = log(m_nu);
  if(v_nu>0){
    m_random_nu = true;
    m_alpha_nu = e_nu*e_nu/v_nu;
    m_beta_nu = m_alpha_nu/e_nu;
  }else
    m_random_nu = false;

  rj_pp_construct();

}

rj_pp::~rj_pp(){

 
    if(m_sample!=NULL){

        destroy_sample();    
        delete [] m_sample;
        m_sample=NULL;
    }
    if(m_current_particle){
        delete m_current_particle;
        m_current_particle = NULL;
    }
   
    delete m_end_of_int_changepoint;
    m_end_of_int_changepoint = NULL;


  if(m_calculate_function){
    delete m_functionofinterest;
  }     
 
  if(m_prop_histogram)
    delete m_prop_histogram;
 
}

void rj_pp::rj_pp_construct(){
  m_spacing_prior = false;
  m_space = 1;
  m_no_neighbouring_empty_intervals = false;
  m_prop_distribution = 'U';
  m_prop_histogram=NULL;
  m_calculate_mean=0;
  m_calculate_function=0;
  m_functionofinterest=NULL;
  m_one_sided_foi = false;
  m_end_of_int_changepoint = new changepoint(m_end_time,0,0,0);
  m_pm->set_data_index(m_end_of_int_changepoint);
  log_m_end_time_minus_start_cps = log(m_end_time-m_start_cps);
}

void rj_pp::proposal_type(const char *proposaltype, void * v){
  
  const char *no = "Normal";
  const char *mix = "Mixture";
  const char *hist = "Histogram";

  if(strcmp(no,proposaltype) == 0){
    m_prop_distribution = 'N';
    m_prop_distribution_mean = ((double*)v)[0];
    m_prop_distribution_sd = ((double*)v)[1];

  }
  else if(strcmp(mix,proposaltype)==0){
    m_prop_distribution = 'M';
    m_prop_distribution_mean = ((double*)v)[0];
    m_prop_distribution_sd = ((double*)v)[1];
  }else if (strcmp(hist,proposaltype)==0 && ((unsigned int*)v)[0]>0){
    m_prop_distribution = 'H';
    m_prop_histogram = new Histogram(m_start_cps,m_end_time,((unsigned int*)v)[0]);
    m_prop_histogram->bin_data(m_pm->get_data());
  }else{
    cerr << "proposal distribution for changepionts is not recognised, using the default uniform distribution" << endl;
  }

}

void rj_pp::initialise_function_of_interest(int grid, bool g, bool prob, bool instant, double delta){
    m_calculate_function=1;   
    m_length_grid=grid;
    m_functionofinterest=new Function_of_Interest(m_length_grid,m_start_time,m_end_time,m_nu,g,prob,m_calculate_mean,0,0,0,instant,delta);

}

changepoint* rj_pp::generate_new_parameter()const {
  double new_value,which;
  which = -1;
  changepoint *cpobj;

  if (m_discrete) {
    bool alreadyexists = true;
    cpobj = new changepoint(0, 0);
    while (alreadyexists) {
      new_value = gsl_rng_uniform_int(r, m_end_time - m_start_cps - 1) + 1 + m_start_cps;
      cpobj -> setchangepoint(new_value);
      alreadyexists = m_current_particle->does_particle_exist(cpobj);
    }
  } else {
    if(m_prop_distribution == 'M'){
      which = gsl_ran_flat( r, 0, 1 );
    }

    if (m_prop_distribution == 'U' || (which<0.5 && which>=0)){
      new_value=gsl_ran_flat( r, m_start_cps, m_end_time );
    }
    else if(m_prop_distribution == 'N' || (which>=0.5)){
      new_value = gsl_ran_gaussian (r, m_prop_distribution_sd);
      new_value += m_prop_distribution_mean;
      double value_theta_right = m_end_time;
      double value_theta_left = m_start_cps;
      double L =  value_theta_right-value_theta_left;
      if (new_value>value_theta_right){
	double delta = new_value - value_theta_right;
	int DL = static_cast<int>(floor(delta/L));
	if (DL % 2 == 0){
	  new_value = value_theta_right - (delta/L-DL)*L;
	}
	else{
	  new_value = value_theta_left + (delta/L - DL)*L;
	}
      }

      if (new_value<value_theta_left){
	double delta = value_theta_left - new_value;
	int DL = static_cast<int>(floor(delta/L));
	if (DL % 2 == 0){           
	  new_value = value_theta_left + ((delta/L)-DL)*L;
	}
	else{
	  new_value = value_theta_right - (delta/L-DL)*L;
	}
      }
    }else if(m_prop_distribution == 'H'){
      new_value = m_prop_histogram->sample_by_differences();
    }else{
      cerr << "can't recognized proposal distribution in RJMCMC_PP.h" << endl;
      exit(1);
    }
    cpobj = new changepoint(new_value,0);
  }
  //    unsigned int new_data_index = m_discrete? get_data_index(new_value,0) : 0; if !m_discrete, delay finding index until neighbouring changepoints are identified.
  
  return(cpobj);
}


double rj_pp::log_likelihood_ratio_birth(changepoint * new_value, int position){

    int k = m_current_particle->get_dim_theta();
    if(m_spacing_prior && k>0){
      if(position<k && m_current_particle->get_theta_component(position)->getchangepoint()-new_value->getchangepoint()<m_space)//check cp to the right
	return -DBL_MAX;
      else if(position>0 && new_value->getchangepoint()-m_current_particle->get_theta_component(position-1)->getchangepoint()<m_space)//check cp to the left
	return -DBL_MAX;
    }

    changepoint *cpobj_right;
    if (position == k){
      cpobj_right = m_end_of_int_changepoint;
      if(!m_conjugate)
	m_pm->propose_new_parameters(m_current_particle,position,0,new_value,cpobj_right);
    }
    else{
      cpobj_right=m_current_particle->get_theta_component(position);
      if(!m_conjugate)
	m_pm->propose_new_parameters(m_current_particle,position,0,new_value,NULL);
    }

    changepoint *cpobj_left = m_current_particle->get_theta_component(position-1);

    /* if there is an empty interval for a discrete time process ignore iteration */
    /*if(m_discrete){
       long long int r=cpobj_right->getdataindex();
       long long int n=new_value->getdataindex();
       long long int l=cpobj_left->getdataindex();
       if((n-l)<=0){
	 return(-1e300);
      }
      if((r-n)<0){
	return(-1e300);
       }
       }*/

    /* if there is an empty interval ignore iteration for cts time process*/
    if(m_no_neighbouring_empty_intervals && !m_discrete){
      unsigned long long int r = cpobj_right->getdataindex();
      unsigned long long int l = cpobj_left->getdataindex();
      if((r-l)<=0)
        return( -1e300);

      if(k>=1){
	double c=new_value->getdataindex();

	if(c-l<=0){
                
	  if(position-1>=0){
                   
	    changepoint *cpobj_double_left=m_current_particle->get_theta_component(position-2);
	    double double_l=cpobj_double_left->getdataindex();
            
	    if((l-double_l<=0))
	      return(-1e300);
	  }
	}else if (r-c<=0){
	  if (position != k){
	    changepoint * cpobj_double_right;
	    if(position+1==k){
	      cpobj_double_right = m_end_of_int_changepoint;
	    }
	    else{
	      cpobj_double_right=m_current_particle->get_theta_component(position+1);
	    }
	    double double_r=cpobj_double_right->getdataindex();
	    if(double_r-c<=0)
	      return(-1e300);
	  }             
	}
      }
    }
    
    //cout << "start " << cpobj_right->getchangepoint() << " " << cpobj_left->getchangepoint() << endl;
    double likelihood_contribution_right = m_pm->log_likelihood_interval(new_value, cpobj_right, cpobj_left);

    new_value -> setlikelihood(likelihood_contribution_right);

 
    if (m_calculate_mean && m_conjugate){
        double mean = m_pm->calculate_mean(new_value, cpobj_right);
        new_value->setmeanvalue(mean);
        new_value->setvarvalue(m_pm->get_var());
    }
   
    double likelihood_contribution_left = m_pm->log_likelihood_interval(cpobj_left,new_value,position>0?m_current_particle->get_theta_component(position-2):NULL);
    
    double old_likelihood = cpobj_left->getlikelihood();
    double log_likelihood_ratio = likelihood_contribution_right + likelihood_contribution_left - old_likelihood;

    return(log_likelihood_ratio);
    
}

double rj_pp::log_likelihood_ratio_death(unsigned int index_theta_delete ){
        
    changepoint *cpobj_left = m_current_particle->get_theta_component(index_theta_delete-1);
    changepoint *cpobj_right;

    if (index_theta_delete == (m_current_particle->get_dim_theta()-1)){
        cpobj_right = m_end_of_int_changepoint;
	if(!m_conjugate)
	  m_pm->propose_new_parameters(m_current_particle,index_theta_delete,1,NULL,m_end_of_int_changepoint);
    }
    else{
        cpobj_right = m_current_particle->get_theta_component(index_theta_delete+1);
	if(!m_conjugate)
	  m_pm->propose_new_parameters(m_current_particle,index_theta_delete,1,NULL,NULL);
    }

    changepoint *cpobj = m_current_particle->get_theta_component(index_theta_delete);
    double likelihood_contribution = m_pm->log_likelihood_interval(cpobj_left,cpobj_right,index_theta_delete>0?m_current_particle->get_theta_component(index_theta_delete-2):NULL);
    double old_likelihood = cpobj_left->getlikelihood()+cpobj->getlikelihood();
    double log_likelihood_ratio = likelihood_contribution - old_likelihood;

    return(log_likelihood_ratio);

}

changepoint* rj_pp::copy_parameter(changepoint *cpcopy) const{
  changepoint * cp = new changepoint(cpcopy);
  return cp;
}

changepoint* rj_pp::move_parameter(unsigned int index_theta_move) const{

    changepoint *cpobj_move = m_current_particle->get_theta_component(index_theta_move);    
    changepoint *cpobj_left = m_current_particle->get_theta_component(index_theta_move-1);
    changepoint *cpobj_right;
    double value_theta_right;

    if (index_theta_move == m_current_particle->get_dim_theta()-1){
        cpobj_right = m_end_of_int_changepoint;
        value_theta_right = m_end_time;
    }
    else{    
        cpobj_right =  m_current_particle->get_theta_component(index_theta_move+1);
        value_theta_right = cpobj_right->getchangepoint();
    }

    double value_theta_move = cpobj_move->getchangepoint();
    double value_theta_left = cpobj_left->getchangepoint();

    if(value_theta_left<m_start_cps){
      value_theta_left=m_start_cps;
    }
    changepoint *cpobj;
    double new_value;
    if (m_discrete) {
      double rightvalue = min((int)value_theta_move + (int) m_move_tolerance, (int)value_theta_right);
      double leftvalue = min(max(0, (int)value_theta_move - (int)m_move_tolerance), (int)value_theta_left);
      if (leftvalue - rightvalue - 2 > 1){
	double newposition = gsl_rng_uniform_int(r, leftvalue - rightvalue - 2);
	new_value = newposition + rightvalue;
	if (new_value == value_theta_move) {
	  new_value += 1;
	}
	cpobj = new changepoint(new_value, (int)new_value);
	return (cpobj);
      } else {
	new_value = value_theta_move;
	cpobj = new changepoint(new_value, cpobj_move->getdataindex());
	return (cpobj);
      }
    }
      
      
    new_value = gsl_ran_flat( r, value_theta_move-m_move_tolerance, value_theta_move+m_move_tolerance );//(2.0*m_move_tolerance)*(double)rand()/(double)RAND_MAX+value_theta_move-m_move_tolerance;
 
   
    double L =  value_theta_right-value_theta_left;

    if (new_value>value_theta_right){

      double delta = new_value - value_theta_right;
      int DL = static_cast<int>(floor(delta/L));

      if (DL % 2 == 0){
	new_value = value_theta_right - (delta/L-DL)*L;
      }
      else{
	new_value = value_theta_left + (delta/L - DL)*L;
      }
    }

    if (new_value<value_theta_left){

      double delta = value_theta_left - new_value;
      int DL = static_cast<int>(floor(delta/L));

      if (DL % 2 == 0){           
	new_value = value_theta_left + ((delta/L)-DL)*L;       
      }
      else{
	new_value = value_theta_right - (delta/L-DL)*L;
      }
    }

    cpobj = new changepoint(new_value,0);
    m_pm->set_data_index(cpobj,0,cpobj_left,cpobj_right);
    return(cpobj);
}

double rj_pp::log_likelihood_ratio_move(changepoint * new_theta, unsigned int k) {
  
    changepoint *cpobj_left = m_current_particle->get_theta_component(k-1);
    changepoint *cpobj_move = m_current_particle->get_theta_component(k);
    changepoint *cpobj_right;
    unsigned int dim = m_current_particle->get_dim_theta();

    if (k == dim-1){
        cpobj_right = m_end_of_int_changepoint;
	if(!m_conjugate)
	  m_pm->propose_new_parameters(m_current_particle,k,2,new_theta,m_end_of_int_changepoint);
    }
    else{    
        cpobj_right =  m_current_particle->get_theta_component(k+1);
	if(!m_conjugate)
	  m_pm->propose_new_parameters(m_current_particle,k,2,new_theta,NULL);
    }

    if(m_spacing_prior && dim>0){
      if(k<dim-1 && cpobj_right->getchangepoint()-new_theta->getchangepoint()<m_space)//check cp to the right
	return -DBL_MAX;
      else if(k>0 && new_theta->getchangepoint()-cpobj_left->getchangepoint()<m_space)//check cp to the left
	return -DBL_MAX;
    }
    
    
    /* if discrete check to make sure no empty intervals*/  
    /*if(m_discrete){
      long long int r=cpobj_right->getdataindex();
      long long int m=new_theta->getdataindex();
      long long int l=cpobj_left->getdataindex();
      long long int move=cpobj_move->getdataindex();

      if((move-m)==0)
	return(-1e300);
 
      if((m-l)<=0)
	return(-1e300);

      if((r-m)<0)
      return(-1e300);

      }*/

    double likelihood_contribution_right = m_pm->log_likelihood_interval(new_theta, cpobj_right, cpobj_left);
    new_theta->setlikelihood(likelihood_contribution_right);

    if (m_calculate_mean && m_conjugate){
        double mean = m_pm->calculate_mean(new_theta, cpobj_right);
        new_theta->setmeanvalue(mean);
        new_theta->setvarvalue(m_pm->get_var());
    }

    double likelihood_contribution_left = m_pm->log_likelihood_interval(cpobj_left,new_theta,k>0?m_current_particle->get_theta_component(k-2):NULL);
    double old_likelihood = cpobj_left->getlikelihood()+cpobj_move->getlikelihood();
    double log_likelihood_ratio = likelihood_contribution_right+likelihood_contribution_left-old_likelihood;

    return(log_likelihood_ratio);
}

double rj_pp::log_likelihood_ratio_move_parameter(int k) {

  changepoint *cpobj_move = m_current_particle->get_theta_component(k);
  changepoint *cpobj_right;

  if (k == (int)(m_current_particle->get_dim_theta()-1)){
    cpobj_right = m_end_of_int_changepoint;
    m_pm->propose_new_parameters(m_current_particle,k,3,NULL,m_end_of_int_changepoint);
  }
  else{    
    cpobj_right =  m_current_particle->get_theta_component(k+1);
    m_pm->propose_new_parameters(m_current_particle,k,3,NULL,NULL);
  }

  double likelihood_contribution_old = cpobj_move->getlikelihood();
  double likelihood_contribution_new = m_pm->log_likelihood_interval(cpobj_move, cpobj_right, m_current_particle->get_theta_component(k-1));
  double log_likelihood_ratio = likelihood_contribution_new-likelihood_contribution_old;

  return(log_likelihood_ratio);
}

double rj_pp::log_prior_ratio_birth(changepoint *newtheta) const{
  double prior=0;

  if (m_discrete) {
    return m_log_nu + log(1-m_nu);
  }
  if(!m_conjugate)
   prior= m_pm->calculate_prior_ratio(m_current_particle,0);

 if(!m_random_nu) {
   if (m_spacing_prior) {
     double distance_to_end=m_end_time - newtheta->getchangepoint();
     prior = m_nu * (distance_to_end<m_space?distance_to_end:m_space);
   }
   return (prior+m_log_nu);
 }

  return (prior+log((m_alpha_nu + m_k)/(m_beta_nu + m_end_time-m_start_cps)));
}

double rj_pp::log_prior_ratio_death(unsigned int thetadelete) const{
  double prior=0;

  if (m_discrete) {
    return - m_log_nu - log(1-m_nu);
  }
  if(!m_conjugate)
    prior=m_pm->calculate_prior_ratio(m_current_particle,1);
  if(!m_random_nu){
    if (m_spacing_prior) {
      double distance_to_end=m_end_time-m_current_particle->get_theta_component(thetadelete)->getchangepoint();
      prior = -m_nu * (distance_to_end<m_space?distance_to_end:m_space);
    }
    return (prior-m_log_nu);
  }
  return (prior+log((m_beta_nu + m_end_time-m_start_cps)/(m_alpha_nu + m_k-1)));
}

double rj_pp::log_prior_ratio_move(changepoint *newposition, unsigned int changepointmove) const{
   double prior=0;
  if(!m_conjugate)
    prior=m_pm->calculate_prior_ratio(m_current_particle,2);  
  if (m_spacing_prior) {
    double old_distance_to_end=m_end_time-m_current_particle->get_theta_component(changepointmove)->getchangepoint();
    double old_space=old_distance_to_end<m_space?old_distance_to_end:m_space;
    double new_distance_to_end=m_end_time-newposition->getchangepoint();
    double new_space=new_distance_to_end<m_space?new_distance_to_end:m_space;
    prior = m_nu * (new_space-old_space);
  }
  return (0+prior);
}

double rj_pp::log_prior_ratio_move_parameter() const{
  return(m_pm->calculate_prior_ratio(m_current_particle,3));
}

double rj_pp::log_proposal_ratio_birth(changepoint* new_value) const{
    int k = m_current_particle->get_dim_theta();

    double value = 0;

    if (m_prop_distribution == 'U'){

      value += - log(k+1) + log_m_end_time_minus_start_cps;

    }else if(m_prop_distribution == 'N'){

        double new_value_theta = new_value->getchangepoint();
        double dummy = log(gsl_ran_gaussian_pdf (new_value_theta-m_prop_distribution_mean, m_prop_distribution_sd));
        value += - log(k+1) - dummy;

    }else if(m_prop_distribution == 'M'){

         double new_value_theta = new_value->getchangepoint();
         double dummy=  -LOG_TWO+log((1.0/(m_end_time-m_start_cps))+gsl_ran_gaussian_pdf(new_value_theta-m_prop_distribution_mean, m_prop_distribution_sd));
	 value = - log(k+1) - dummy;

    }else if(m_prop_distribution == 'H'){

      double new_value_theta = new_value->getchangepoint();
      value+=-log(k+1)-m_prop_histogram->sampling_by_differences_log_density(new_value_theta);

    }else{
      cerr << "can't recognise proposal distribution in log_prior_ratio_move in RJMCMC_PP.h" << endl;
      exit(1);
    }

    if(!m_conjugate)
      value+=m_pm->proposal_ratio(m_current_particle,0);  

    return(value);
}     

double rj_pp::log_proposal_ratio_death(int index) const{

    int k = m_current_particle->get_dim_theta();

    double value=0;

    if (m_prop_distribution == 'U'){

           value += - log_m_end_time_minus_start_cps + log(k);
       
    }else if(m_prop_distribution == 'N'){

        changepoint * cp = m_current_particle->get_theta_component(index);
        double value_theta_delete = cp->getchangepoint();
        double dummy = log(gsl_ran_gaussian_pdf (value_theta_delete-m_prop_distribution_mean, m_prop_distribution_sd));
        value += dummy + log(k);
       
    }else if(m_prop_distribution == 'M'){
        
        changepoint * cp = m_current_particle->get_theta_component(index);
        double value_theta_delete = cp->getchangepoint();
	double dummy=  -LOG_TWO+log((1.0/(m_end_time-m_start_time))+ gsl_ran_gaussian_pdf(value_theta_delete-m_prop_distribution_mean, m_prop_distribution_sd));
	value += dummy + log(k);
  
    }else if(m_prop_distribution == 'H'){

      changepoint * cp = m_current_particle->get_theta_component(index);
      double value_theta_delete = cp->getchangepoint();
      value+= log(k)+m_prop_histogram->sampling_by_differences_log_density(value_theta_delete);
    }else{
      cerr << "can't recognise proposal distribution in log_prior_ratio_death in RJMCMC_PP.h" << endl;
      exit(1);
    }

    if(!m_conjugate)
      value+=m_pm->proposal_ratio(m_current_particle,1);  

    return(value);
}

double rj_pp::log_proposal_ratio_move() const{
  double value=0;

 if(!m_conjugate)
      value+=m_pm->proposal_ratio(m_current_particle,2);  

 return(value);
}

double rj_pp::log_proposal_ratio_move_parameter() const{

  return(m_pm->proposal_ratio(m_current_particle,3));

}

void rj_pp::update_likelihood_birth(changepoint* new_theta ,int position,double like_ratio){

    changepoint * cpobj_left =  m_current_particle ->get_theta_component(position-1);
    double like_new_left = like_ratio + cpobj_left->getlikelihood() -new_theta->getlikelihood();
    cpobj_left->setlikelihood(like_new_left);

    if (m_calculate_mean && m_conjugate){
        double mean = m_pm->calculate_mean(cpobj_left, new_theta);
        cpobj_left->setmeanvalue(mean);
        cpobj_left->setvarvalue(m_pm->get_var());
    }
}

void rj_pp::update_likelihood_death(unsigned int index, double like_ratio){

    changepoint * cpobj_delete = m_current_particle ->get_theta_component(index);
    changepoint *cpobj_left = m_current_particle -> get_theta_component(index-1);
    
    double like_new = like_ratio +  cpobj_delete->getlikelihood() + cpobj_left->getlikelihood();

    cpobj_left->setlikelihood(like_new);

    if (m_calculate_mean && m_conjugate){

        changepoint * cpobj_right;
        if(index == m_current_particle->get_dim_theta()-1){
            cpobj_right = m_end_of_int_changepoint;
        }
        else{            
            cpobj_right = m_current_particle->get_theta_component(index+1);
        }

        double mean = m_pm->calculate_mean(cpobj_left, cpobj_right);
        cpobj_left->setmeanvalue(mean);
        cpobj_left->setvarvalue(m_pm->get_var());
    }
}

void rj_pp::update_likelihood_move(changepoint * new_theta, int position, double like_ratio){

    changepoint *cpobj_move = m_current_particle->get_theta_component(position);
    changepoint *cpobj_left = m_current_particle->get_theta_component(position-1);

    double like_new = like_ratio - new_theta->getlikelihood() + cpobj_move->getlikelihood() + cpobj_left->getlikelihood();

    cpobj_left->setlikelihood(like_new);

    if (m_calculate_mean && m_conjugate){
        double mean = m_pm->calculate_mean(cpobj_left, new_theta);
        cpobj_left->setmeanvalue(mean);
        cpobj_left->setvarvalue(m_pm->get_var());
    }
}

void rj_pp::update_likelihood_move_parameter(int position, double like_ratio){

    changepoint *cpobj_move = m_current_particle->get_theta_component(position);
    double like_new = like_ratio + cpobj_move->getlikelihood();
    cpobj_move->setlikelihood(like_new);
}



void rj_pp::initiate_sample(Particle<changepoint>* ptr2particle){
 
    if (ptr2particle==NULL){    
      
      changepoint *cpobj = new changepoint(m_start_time,0);
      m_pm->set_data_index(cpobj,0);
      m_current_particle = new Particle<changepoint>(0,0,cpobj);
      changepoint *cpobj_temp =  m_end_of_int_changepoint;

      if(!m_conjugate){
	m_pm->propose_new_parameters(m_current_particle,-1,3,NULL,m_end_of_int_changepoint);
      }

      double likelihood = m_pm->log_likelihood_interval(cpobj,cpobj_temp);
      double prior = -m_nu*(m_end_time-m_start_time);
      cpobj->setlikelihood(likelihood);
      m_log_likelihood = likelihood;
      m_current_particle -> set_log_posterior(likelihood+prior);

      if (m_calculate_mean && m_conjugate){

        double mean = m_pm->calculate_mean(cpobj, cpobj_temp);
        cpobj->setmeanvalue(mean);
        cpobj->setvarvalue(m_pm->get_var());
      }
    }
    else{
      m_current_particle = ptr2particle;
    }
}

void rj_pp::calculate_function_of_interest(){
  double sum=0;
  probability_model * pm = NULL;
  if(!m_conjugate||m_pm->windowed_model())
    pm = m_pm;

  if(m_storing_sample)
    m_functionofinterest->calculate_function(m_start_time,m_end_time,m_sample,m_size_of_sample,NULL,sum,0,0,0,pm);
  else{
    double weight = m_current_importance_weight;
    m_functionofinterest->calculate_function(m_start_time,m_end_time,&m_current_particle,1,&weight,sum,0,0,0,pm);
  }
}

void rj_pp::update_primary_function_of_interest(){
  double mean_i=0, var_i=0, grid_i=0;
  int position_i = -1, last_position = -1;
  bool interior_grid = m_length_grid < 20;
  double cp_left_position = m_start_cps;//m_current_particle->get_theta_component(position_i)->getchangepoint();
  double cp_right_position = (position_i < (int)m_k-1)?m_current_particle->get_theta_component(position_i+1)->getchangepoint():m_end_time;
  for(unsigned int i=0; i<m_length_grid; i++){
    grid_i = m_start_cps + (m_end_time - m_start_cps)*(interior_grid?(i+1)/(double)(m_length_grid+1):i/(double)(m_length_grid-1));
    while( position_i < (int)m_k-1 && cp_right_position< grid_i ){
      position_i++;
      cp_left_position = cp_right_position;
      cp_right_position = (position_i < (int)m_k-1) ? m_current_particle->get_theta_component(position_i+1)->getchangepoint():m_end_time;
    }

    if(m_calculate_mean){
      if(position_i != last_position || i== 0){
	last_position = position_i;

	mean_i = m_current_particle->get_theta_component(position_i)->getmeanvalue();
	var_i = m_current_particle->get_theta_component(position_i)->getvarvalue();
      }
    }
    if(!m_calculate_mean){
      double dist_left = grid_i - cp_left_position;
      double dist_right = cp_right_position - grid_i;
      if(dist_left<dist_right||m_one_sided_foi)
	mean_i = dist_left;
      else
	mean_i = dist_right;
    }
    
    double mean_fn_i = mean_i * m_pm->get_mean_function(grid_i-cp_left_position);
    m_mean_function_of_interest[i] += mean_fn_i;
    if(m_mean_sq_function_of_interest)
      m_mean_sq_function_of_interest[i] += mean_fn_i * mean_fn_i;
    if(m_var_function_of_interest)
      m_var_function_of_interest[i] += var_i;
  }
  
}

void rj_pp::calculate_means( Particle< changepoint> * p){
  if( !p )
    p = get_MAP_dimension_MAP();
  int k = p->get_dim_theta();
  changepoint *cpobj_left = p->get_theta_component(-1);
  changepoint *cpobj_right = NULL;
  for( int position = 0; position <= k; position++ ){
    if(position == k)
      cpobj_right = m_end_of_int_changepoint;
    else
      cpobj_right = p->get_theta_component(position);
    double mean = m_pm->calculate_mean(cpobj_left, cpobj_right);
    cpobj_left->setmeanvalue(mean);
    cpobj_left->setvarvalue(m_pm->get_var());
    cpobj_left = cpobj_right;
  }
}

bool rj_pp::draw_means_from_posterior(bool test_monotonicity){
  changepoint *cpobj_left = m_current_particle->get_theta_component(-1);
  double cp_left_value=cpobj_left->getchangepoint();
  changepoint *cpobj_right = NULL;
  changepoint *cpobj_left2 = NULL;
  bool monotonic = test_monotonicity;
  double last_right_mean = -DBL_MAX;
  for( unsigned int position = 0; position <= m_k; position++ ){
    if(position == m_k)
      cpobj_right = m_end_of_int_changepoint;
    else
      cpobj_right = m_current_particle->get_theta_component(position);
    double cp_right_value=cpobj_right->getchangepoint();
    double mean = m_pm->draw_mean_from_posterior(cpobj_left, cpobj_right,cpobj_left2);
    cpobj_left->setmeanvalue(mean);
    cpobj_left->setvarvalue(m_pm->get_var());
    if(monotonic){
      if(mean<last_right_mean)
	monotonic=false;
      else
	last_right_mean=mean*m_pm->get_mean_function(cp_right_value-cp_left_value);
    }
    cpobj_left2 = cpobj_left;
    cpobj_left = cpobj_right;
    cp_left_value=cp_right_value;
  }
  return monotonic;
}

Step_Function* rj_pp::create_step_function_of_means( Particle< changepoint> * p, bool normalise){
  if( !p )
    p = get_MAP_dimension_MAP();
  unsigned int num_knots = p->get_dim_theta() + 1;
  vector<double> knots, heights;
  for( int i=0; i<(int)num_knots; i++){
    knots.push_back((i==0)?m_start_cps:p->get_theta_component(i-1)->getchangepoint());
    heights.push_back(p->get_theta_component(i-1)->getmeanvalue());
  }
  Step_Function* s = new Step_Function(&knots,&heights,m_end_time,true);
  if(normalise)
    s->standardise_step_function();
  return s;
}

void rj_pp::write_means_to_file( const char* output_filename, Particle< changepoint> * p, bool normalise){
  ofstream OutputStream(output_filename, ios::out);
  write_means_to_ostream(OutputStream,p,normalise);
  OutputStream.close();
}

void rj_pp::write_means_to_ostream( ostream& OutputStream, Particle< changepoint> * p, bool normalise){
  Step_Function* s = create_step_function_of_means( p, normalise );
  OutputStream << s;
  delete s;
}

/*
{
  if( !p )
    p = get_MAP_dimension_MAP();//m_current_particle;
  OutputStream<<setiosflags(ios::fixed);
  OutputStream.precision(16);
  double current_cp = 0, last_cp = m_start_cps, mean_value = 0;
  double scaling_factor = 1;
  if(normalise){
    scaling_factor=0;
    for( int i=0; i<=(int)(p->get_dim_theta()); i++){
      current_cp = i<(int)(p->get_dim_theta())? (p->get_theta_component(i)->getchangepoint()): m_end_time;
      mean_value = p->get_theta_component(i-1)->getmeanvalue();
      scaling_factor += (current_cp-last_cp)*mean_value;
      last_cp = current_cp;
    }
  }
  double integral_value = 0;
  mean_value = 0;
  last_cp = m_start_cps;
  for( int i=-1; i<(int)(p->get_dim_theta()); i++){
    current_cp = (i<0)?m_start_cps:p->get_theta_component(i)->getchangepoint();
    if(i>=0)
      integral_value += (current_cp-last_cp)*mean_value;
    mean_value = p->get_theta_component(i)->getmeanvalue()/scaling_factor;
    OutputStream<< current_cp << "\t"<< mean_value << "\t"<< integral_value <<endl;
    last_cp = current_cp;
  }
}
*/

void rj_pp::write_changepoints_to_file( ostream& CPStream, Particle< changepoint> * p){
  double current_cp;
  
if( !p ){
    p = get_MAP_dimension_MAP();//m_current_particle;
  }

  int size = p->get_dim_theta();
  
  CPStream<<size<<" ";
  for( int i=-1; i<size; i++){
    current_cp = p->get_theta_component(i)->getchangepoint();
    CPStream << current_cp << " ";
  }
  CPStream<<endl;
}

void rj_pp::write_changepoints_and_likelihood_to_file( ostream& CPStream, Particle< changepoint> * p){
  double current_cp;
  double log_lik;
  if( !p ){
    p = get_MAP_dimension_MAP();
  }

  int size = p->get_dim_theta();
  
  CPStream<<size<<" ";
  for( int i=-1; i<size; i++){
    current_cp = p->get_theta_component(i)->getchangepoint();
    log_lik = p->get_theta_component(i)->getlikelihood();
    CPStream << current_cp << " " <<log_lik<<" ";
  }
  CPStream<<endl;
}

void rj_pp::update_new_parameter_with_position( changepoint* cp, unsigned int position ){

  if(position == m_current_particle->get_dim_theta()){
    m_pm->set_data_index(cp,0,m_current_particle->get_theta_component(position-1));
  }else{
    m_pm->set_data_index(cp,0,m_current_particle->get_theta_component(position-1),m_current_particle->get_theta_component(position));
  }
}

void rj_pp::locally_perfect_MAP(Particle<changepoint>* p, bool smoothing){
  if(!p)
    p = get_MAP_dimension_MAP();
  Data<double> * X = m_pm->get_data();
  unsigned int size = p->get_dim_theta();
  double epsilon = 1e-8*(m_end_time-m_start_time);//X->get_range_width();
  for(unsigned int i = 0; i< size;i++){
    changepoint* cpobj_i = p->get_theta_component(i);
    if(m_discrete){
      cpobj_i->setchangepoint(ceil(cpobj_i->getchangepoint()));
      break;
    }
    unsigned long long int index_i = cpobj_i->getdataindex();
    changepoint* cpobj_left = p->get_theta_component(i-1);
    changepoint* cpobj_right = i<size-1 ? p->get_theta_component(i+1) : m_end_of_int_changepoint;
    double old_likelihood = cpobj_left->getlikelihood()+cpobj_i->getlikelihood();
    bool moved_left = false;
    if(index_i > 0){
      double left_value = (*X)[0][index_i-1]+epsilon;// : m_pm->get_start();
      if( i==0 || left_value > p->get_theta_component(i-1)->getchangepoint()){
	changepoint* proposed_cpobj = new changepoint(left_value,index_i);
	double likelihood_contribution_left = m_pm->log_likelihood_interval(cpobj_left,proposed_cpobj,i>0?p->get_theta_component(i-2):NULL);
	double likelihood_contribution_right = m_pm->log_likelihood_interval(proposed_cpobj, cpobj_right, cpobj_left);
	double log_posterior_ratio = likelihood_contribution_left+likelihood_contribution_right-old_likelihood;
	if(log_posterior_ratio>0){
	  proposed_cpobj->setlikelihood(likelihood_contribution_right);
	  p->change_component(proposed_cpobj,i);
	  cpobj_left->setlikelihood(likelihood_contribution_left);
	  double left_mean = m_pm->calculate_mean(cpobj_left, proposed_cpobj);
	  cpobj_left->setmeanvalue(left_mean);
	  cpobj_left->setvarvalue(m_pm->get_var());
	  double mean = m_pm->calculate_mean(proposed_cpobj, cpobj_right);
	  proposed_cpobj->setmeanvalue(mean);
	  proposed_cpobj->setvarvalue(m_pm->get_var());
	  moved_left = true;
	  if(smoothing){
	    double tau1 = (*X)[0][index_i-1];
	    double tau2 = index_i < X->get_cols() ? (*X)[0][index_i]:m_end_time;
	    //	    double tau = tau1+(tau2-tau1)*(mean)/(left_mean+mean);
	    double tau = left_mean > mean ? (tau2*exp((mean-left_mean)*(tau2-tau1))-tau1)/(exp((mean-left_mean)*(tau2-tau1))-1) : (tau2-tau1*exp((mean-left_mean)*(tau1-tau2)))/(1-exp((mean-left_mean)*(tau1-tau2)));
	    tau -= 1/(mean-left_mean);
	    proposed_cpobj->setchangepoint(tau);
	  }
	}
	else
	  delete proposed_cpobj;
      }
    }
    if(index_i < X->get_cols() && !moved_left){
      double right_value = (*X)[0][index_i];// - epsilon;// : m_pm->get_end();
      if( i==size-1 || right_value < p->get_theta_component(i+1)->getchangepoint()){
	changepoint* proposed_cpobj = new changepoint(right_value,index_i);
	double likelihood_contribution_left = m_pm->log_likelihood_interval(cpobj_left,proposed_cpobj,i>0?p->get_theta_component(i-2):NULL);
	double likelihood_contribution_right = m_pm->log_likelihood_interval(proposed_cpobj, cpobj_right, cpobj_left);
	proposed_cpobj->setlikelihood(likelihood_contribution_right);
	double log_posterior_ratio = likelihood_contribution_left+likelihood_contribution_right-old_likelihood;
	if(log_posterior_ratio>0){
	  proposed_cpobj->setlikelihood(likelihood_contribution_right);
	  p->change_component(proposed_cpobj,i);
	  cpobj_left->setlikelihood(likelihood_contribution_left);
	  double left_mean = m_pm->calculate_mean(cpobj_left, proposed_cpobj);
	  cpobj_left->setmeanvalue(left_mean);
	  cpobj_left->setvarvalue(m_pm->get_var());
	  double mean = m_pm->calculate_mean(proposed_cpobj, cpobj_right);
	  proposed_cpobj->setmeanvalue(mean);
	  proposed_cpobj->setvarvalue(m_pm->get_var());
	  if(smoothing){
	    double tau1 = index_i > 0 ? (*X)[0][index_i-1] : m_start_time;
	    double tau2 = (*X)[0][index_i];
	    //	    double tau = tau1+(tau2-tau1)*(mean)/(left_mean+mean);
	    double tau = left_mean > mean ? (tau2*exp((mean-left_mean)*(tau2-tau1))-tau1)/(exp((mean-left_mean)*(tau2-tau1))-1) : (tau2-tau1*exp((mean-left_mean)*(tau1-tau2)))/(1-exp((mean-left_mean)*(tau1-tau2)));
	    tau -= 1/(mean-left_mean);
	    proposed_cpobj->setchangepoint(tau);
	  }
	}
	else
	  delete proposed_cpobj;
      }
    }
  }
}


void rj_pp::calculate_importance_weight(){

  m_current_log_importance_weight=0;


  m_pm->set_prior_parameters(NULL,m_current_particle->get_theta_component(-1));
  double prior_parameter_1 = m_pm->get_alpha();
  double prior_parameter_2 = m_pm->get_beta();
  double intep=(m_current_particle->get_theta_component(-1))->getmeanvalue(); 	 
  m_current_log_importance_weight+=log(gsl_ran_gamma_pdf (intep, 4.5, (double)1.0/1.5));
  m_current_log_importance_weight-=log(gsl_ran_gamma_pdf (intep, prior_parameter_1, (double)1.0/prior_parameter_2));

//m_intensity[k]/=gsl_ran_gamma_pdf (inte, 0.1, (double)1.0/0.1);
  for(unsigned int i=0; i<m_current_particle->get_dim_theta(); i++){
    m_pm->set_prior_parameters(m_current_particle->get_theta_component(i-1),m_current_particle->get_theta_component(i));
    prior_parameter_1 = m_pm->get_alpha();
    prior_parameter_2 = m_pm->get_beta();
    double inte=(m_current_particle->get_theta_component(i))->getmeanvalue();
    intep=(m_current_particle->get_theta_component(i-1))->getmeanvalue();
      m_current_log_importance_weight+=log(gsl_ran_gamma_pdf (inte, intep*intep/5, (double)5/intep));
      m_current_log_importance_weight-=log(gsl_ran_gamma_pdf (inte, prior_parameter_1, (double)1.0/prior_parameter_2));
  }
}
