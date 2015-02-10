#ifndef FUNCTION_OF_INTEREST_H
#define FUNCTION_OF_INTEREST_H

#include "particle.h"
#include "changepoint.h"
#include "probability_model.h"

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
  double get_variance_distance(){return m_variance_distance;}
  double calculate_variance_prob(int);
  void reset_prob();
  void set_start(int st){m_start_of_sample = st;}
  void write_mean_to_file(const string output_filename = "intensity.txt");
  void set_importance_sampling(){m_coal_importance_sampling = 1;}

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
   long double m_average_distance_squared;
   double m_variance_distance;
   int m_start_of_sample;
   int m_sample_size;
   
   //do importance sampling to compare with Del Moral as in the smc paper
   bool m_coal_importance_sampling;

};

Function_of_Interest::Function_of_Interest(int grid,double start, double end, double prior_term, bool g, bool prob, bool intensity, bool online, bool sequential,int intervals,bool instant, double delta )
:m_grid(grid),m_start(start),m_end(end),m_prior(prior_term),m_calculate_g(g),m_calculate_prob(prob),m_calculate_intensity(intensity),m_online(online),m_update(sequential),m_instantaneous(instant),m_delta(delta)
{

  m_coal_importance_sampling=0;
   m_prior_expectation_function_of_interest=NULL;
   m_prior_sd_function_of_interest=NULL;
   m_exp_last_changepoint_sequential=NULL;
   m_prob_last_changepoint_sequential=NULL;
   m_exp_last_changepoint=NULL;
   m_variance_exp_last_changepoint=NULL;
   m_prob_function_of_interest=NULL;
   m_intensity=NULL;


  m_start_of_sample=0;

  if (m_online){
    m_fixed=0;
  }
  else{
    m_fixed=1;
  }
  
  m_grid_points =(end-start)/static_cast<double>(grid);
   
  m_cont=0;
 
  if (m_calculate_g){
    if(m_online || m_fixed){
      m_exp_last_changepoint = new long double[m_grid];
      m_variance_exp_last_changepoint = new long double[m_grid];
     
    }
        
    if(m_update){
      m_exp_last_changepoint_sequential = new double *[static_cast<int>(intervals)];
      m_exp_last_changepoint_sequential[0] = new double[m_grid*(static_cast<int>(intervals))];
      for (int i=1; i<intervals; i++){
	m_exp_last_changepoint_sequential[i]=m_exp_last_changepoint_sequential[i-1]+m_grid;
      }
    }
  }


       
  if (m_calculate_prob){
    if(m_online || m_fixed){
      m_prob_function_of_interest = new long double[m_grid];
	
    }
    if(m_update){
      m_prob_last_changepoint_sequential = new double *[static_cast<int>(intervals)];
      m_prob_last_changepoint_sequential[0] = new double[m_grid*(static_cast<int>(intervals))];
      for (int i=1; i<intervals; i++){
	m_prob_last_changepoint_sequential[i]=m_prob_last_changepoint_sequential[i-1]+m_grid;
      }
    }
  }


  for(int i=0; i<m_grid; i++){
    if(m_online || m_fixed){
      if(m_calculate_g){  
	m_exp_last_changepoint[i]=0;
	m_variance_exp_last_changepoint[i]=0;
      }
      if(m_calculate_prob)
	m_prob_function_of_interest[i]=0;
    }
    if(m_update){
      for(int j=0; j<intervals; j++){
	if (m_calculate_g)
	  m_exp_last_changepoint_sequential[j][i]=0;
	if(m_calculate_prob)
	  m_prob_last_changepoint_sequential[j][i]=0;
      }
    }
  }
    

  if(m_calculate_intensity){

    m_intensity = new double[m_grid];
    for(int i=0; i<m_grid; i++)
      m_intensity[i]=0;
  }
    
  double dummy_variable;
  bool test;

  m_prior_expectation_function_of_interest = new double[m_grid];
  m_prior_sd_function_of_interest = new double[m_grid];

  for(int i=0; i<m_grid; i++){

    if(m_cont){
      (m_cont && (m_grid_points*(i+1)<(m_end-m_delta) && m_grid_points*(i+1)>m_delta) ? test=1 : test=0);
    }else{
      (m_grid_points*(i+1)>m_delta ? test=1 : test=0);
    }
    if(m_instantaneous && test){
      dummy_variable=m_prior*m_delta;
    }
    else{
      if(m_cont && m_grid_points*(i+1)>((m_end-m_start)/2)){
	dummy_variable=m_prior*(m_end-m_grid_points*(i+1));
      }else{
	dummy_variable=m_prior*m_grid_points*(i+1);
      }
    }

  
    if(m_cont){
      m_prior_expectation_function_of_interest[i]=(1.0/(2*m_prior))*(1-exp(-2*dummy_variable));
      m_prior_sd_function_of_interest[i]=((1.0/(4*m_prior*m_prior))*(2.0)*(1-exp(-2*dummy_variable)*(1+2*dummy_variable)))-m_prior_expectation_function_of_interest[i]*m_prior_expectation_function_of_interest[i];
    }else{
      m_prior_expectation_function_of_interest[i]=(1.0/m_prior)*(1-exp(-dummy_variable));
      m_prior_sd_function_of_interest[i]=((1.0/(m_prior*m_prior))*(2.0)*(1-exp(-dummy_variable)*(1+dummy_variable)))-m_prior_expectation_function_of_interest[i]*m_prior_expectation_function_of_interest[i];
    }
      
    m_prior_sd_function_of_interest[i]=sqrt(m_prior_sd_function_of_interest[i]);

    
    //squared function of interest
    /*		m_prior_expectation_function_of_interest[i]=(1.0/(prior*prior))*(2.0)*exp(-dummy_variable)*(exp(dummy_variable)-1-dummy_variable);


		m_prior_sd_function_of_interest[i]=((24.0/(prior*prior*prior*prior))*exp(-dummy_variable)*(exp(dummy_variable)-(1+dummy_variable+((dummy_variable*dummy_variable)/2.0)+((dummy_variable*dummy_variable*dummy_variable)/6.0))))-m_prior_expectation_function_of_interest[i]*m_prior_expectation_function_of_interest[i];


		m_prior_sd_function_of_interest[i]=sqrt(m_prior_sd_function_of_interest[i]);*/

  }
}


Function_of_Interest::~Function_of_Interest(){

  if (m_calculate_g){
    if(m_online || m_fixed){
      delete [] m_exp_last_changepoint;
      delete [] m_variance_exp_last_changepoint;
    }
      
    if(m_update){

      delete [] m_exp_last_changepoint_sequential[0];
      delete [] m_exp_last_changepoint_sequential;
    }
  }

       
  if (m_calculate_prob){
    if(m_online || m_fixed){
      delete [] m_prob_function_of_interest;

    }
    if(m_update){
      delete [] m_prob_last_changepoint_sequential[0];
      delete [] m_prob_last_changepoint_sequential;
    }
  }


  if(m_calculate_intensity){
    delete [] m_intensity;
  }
 
  delete [] m_prior_expectation_function_of_interest;
  delete [] m_prior_sd_function_of_interest;
}

void Function_of_Interest::reset_prob(){
 
  for(int i=0; i<m_grid; i++){
 
    if(m_prob_function_of_interest){
      m_prob_function_of_interest[i]=0;}
    if(m_exp_last_changepoint)
      m_exp_last_changepoint[i]=0;
    if(m_intensity)
      m_intensity[i]=0;

  }

  /*for(int i=0; i<m_grid; i++){
  m_prob_function_of_interest[i]*=m_sample_size;
  }*/
}

double Function_of_Interest::calculate_variance_prob(int sample_size){

  double variance=0;
  double variance_future=0;
  
  for(int i=0; i<m_grid; i++){
    variance+= m_prob_function_of_interest[i]*(1-m_prob_function_of_interest[i])/(sample_size-1);
    variance_future+=(sample_size)*(m_prob_function_of_interest[i]*(1-m_prob_function_of_interest[i]))/(double)((sample_size-1)*(sample_size+10));
  }
  variance/=(int)m_grid;
  variance_future/=(int)m_grid;
  
  return(variance);//-variance_future);
}
  
   
void Function_of_Interest::calculate_function(double interval_begin, double interval_end,Particle<changepoint> ** sample,long long int sample_size, double * weights, double sum_weights, double sum_weights_squared, int iters,bool normalise, probability_model * pm){
  double begin;

  bool create_weights=0;
  if(!weights)
    create_weights=1;
  
  m_sample_size=sample_size;

  if(m_update || m_fixed)
    begin=m_start; 
  else
    begin=interval_begin;


  if(m_coal_importance_sampling && !m_fixed){
    Particle<changepoint> * temp;
    long double temp_weight;
    sum_weights=0;
    for(int j=0; j<sample_size; j++){
      temp=sample[j];
      temp_weight=0;
      double prior_parameter_1 = pm->get_alpha();
      double prior_parameter_2 = pm->get_beta();
      double intep=(temp->get_theta_component(-1))->getmeanvalue(); 	 
      temp_weight+=log(gsl_ran_gamma_pdf (intep, 4.5, (double)1.0/1.5));
      temp_weight-=log(gsl_ran_gamma_pdf (intep, prior_parameter_1, (double)1.0/prior_parameter_2));

      for(unsigned int i=0; i<temp->get_dim_theta(); i++){
	double inte=(temp->get_theta_component(i))->getmeanvalue();
	intep=(temp->get_theta_component(i-1))->getmeanvalue();
	temp_weight+=log(gsl_ran_gamma_pdf (inte, intep*intep/5, (double)5/intep));
	temp_weight-=log(gsl_ran_gamma_pdf (inte, prior_parameter_1, (double)1.0/prior_parameter_2));
      }
      weights[j]*=exp(temp_weight);
      sum_weights+=weights[j];
    }
  }

  
  if (create_weights){
    weights=new double[sample_size];
    for(int i=0; i< sample_size; i++)
      weights[i]=1;
    sum_weights=sample_size;
    sum_weights_squared=sample_size;
  }
  
  if(m_online){
    m_average_distance=0;
    m_average_distance_squared=0;
    m_variance_distance=0;
  }

  long double g ,g1=0;
  int loc_index_1,loc_index_2,start_index,ind,dim;
  changepoint * cpobj;
  changepoint * cpobj1=NULL;
  changepoint * end_changepoint=NULL;
  if(m_cont){
    end_changepoint = new changepoint(interval_end,0,0,0);
  }
  
  int fb=static_cast<int>(floor((10000*(double)interval_begin)/(10000*(double)m_grid_points)));
  int fe =static_cast<int>(floor((10000*(double)interval_end)/(10000*(double)m_grid_points)));

  
  for (int i=m_start_of_sample; i<sample_size; i++){

    loc_index_1=static_cast<int>(floor((10000*begin)/(10000*m_grid_points)));
   
    dim = sample[i]->get_dim_theta();
    
    ind=0;
    start_index=-1;   	


    while ( dim>0 && ind == 0 ){
      for (int j=0; j<dim; j++){
	cpobj = sample[i]->get_theta_component(j);

	if(cpobj->getchangepoint()>(ceil((10000*begin)/(10000*m_grid_points))*m_grid_points)){
	  ind=1;
	  start_index=j-1;
	  break;
	}
	if(j==(dim-1)){
	  start_index=dim-1;
	  ind=1;
	}
      }
    }
   
    for (int j=start_index; j<dim; j++){
      if(j==(dim-1)){                
	loc_index_2=fe;
      }
      else{
	cpobj = sample[i]->get_theta_component(j+1);
	loc_index_2 = static_cast<int>(floor((10000*((double)cpobj->getchangepoint()))/(10000*(double)m_grid_points)));
      }

      
      cpobj=sample[i]->get_theta_component(j);

      if(m_cont){
	if(j<(dim-1)){
	  cpobj1=sample[i]->get_theta_component(j+1);
	}
	else{
	  cpobj1=end_changepoint;
	}
      }


      for (int k=loc_index_1; k<loc_index_2; k++){
	
	g= (m_grid_points*(k+1) - cpobj->getchangepoint());
	
	if(m_cont){
	  g1=cpobj1->getchangepoint()-(m_grid_points*(k+1));
	}

	
	if(m_online && k==(fe-1)){
	  m_average_distance+=g*(weights[i]/sum_weights);
	  //m_average_distance_squared+=g*g*weights[i];
	}

	if(m_calculate_intensity && k>=fb){
	  double e;
	  pm?e = pm->get_mean_function(g):e=1;
	  m_intensity[k] += cpobj->getmeanvalue()*e*weights[i];
   	}
	
	if(m_instantaneous){
	  if(g>m_delta)
	    g=m_delta;
	  if(m_cont){
	    if(g1>m_delta)
	      g1=m_delta;
	  }
	}

	if(m_cont){
	  if(g1<=g)
	    g=g1;
	}
	
	if(m_calculate_g){
	  if(m_update)
	    m_exp_last_changepoint_sequential[iters][k]+=g*weights[i];
	  
	  if(k>=fb){
	    m_exp_last_changepoint[k] += g*weights[i];
	    m_variance_exp_last_changepoint[k]+=g*g*(weights[i]/(double)sum_weights);
	  }
	}
	    

	if (m_calculate_prob){
	  if(m_update && (g-m_prior_expectation_function_of_interest[k])<0){
	    m_prob_last_changepoint_sequential[iters][k]+=weights[i];
	  }

	  if (k>=fb && (g-m_prior_expectation_function_of_interest[k])<0){
	    m_prob_function_of_interest[k]+=weights[i];
	  }
	}    

	loc_index_1=loc_index_2;

      }
    }
  }

 
  loc_index_1=(int)floor((10000*begin)/(10000*m_grid_points));
      
  loc_index_2=fe;


  /*if(m_online){
  // m_average_distance/=sum_weights;
  m_average_distance_squared/=sum_weights;
  m_variance_distance=m_average_distance_squared-m_average_distance*m_average_distance;
  }*/
  if(m_online){
    if(m_average_distance>interval_end){
      m_average_distance=interval_end;
    }
  }

  if(normalise){
    normalise_function(sum_weights,loc_index_1,loc_index_2,fb,iters);
  }
 
  if(create_weights){
    delete [] weights; 
  } 

  if(m_cont){
    delete end_changepoint;
  }
 
}

void Function_of_Interest::normalise_function(double sum_weights, int loc_index_1, int loc_index_2, int fb,int iters){


  for (int i=loc_index_1; i<loc_index_2; i++){
      
    if(m_update){
      if(m_calculate_g){
	m_exp_last_changepoint_sequential[iters][i] /= sum_weights;
	// m_exp_last_changepoint_sequential[iters][i] -= m_prior_expectation_function_of_interest[i];
	//m_exp_last_changepoint_sequential[iters][i] /= m_prior_sd_function_of_interest[i]; 
      }
      if(m_calculate_prob){
	m_prob_last_changepoint_sequential[iters][i] /= sum_weights;
      }
    }
	  
    if(i>=fb){

      if(m_calculate_g){
	m_exp_last_changepoint[i]/=sum_weights;
	//m_variance_exp_last_changepoint[i]*=sum_weights;
	m_variance_exp_last_changepoint[i]-=m_exp_last_changepoint[i]*m_exp_last_changepoint[i];
	//m_variance_exp_last_changepoint[i]/=(double)(sum_weights*sum_weights-sum_weights_squared);
	//  m_variance_exp_last_changepoint[i]/=(m_prior_sd_function_of_interest[i]*m_prior_sd_function_of_interest[i]);
	m_exp_last_changepoint[i] -= m_prior_expectation_function_of_interest[i];
	m_exp_last_changepoint[i]/= m_prior_sd_function_of_interest[i];
      }
      if(m_calculate_prob){
	m_prob_function_of_interest[i] /=  sum_weights;
      }
      if(m_calculate_intensity){
	m_intensity[i] = m_intensity[i]/sum_weights;
      }	  
    }
  }
}


void Function_of_Interest::write_mean_to_file(const string output_filename){
  if(!m_intensity){
    cerr << "function_of_interest.h: intensity has not been calculated" << endl;
    return;
  }

  ofstream OutputStream(output_filename.c_str(), ios::out);
  OutputStream << setiosflags(ios::fixed);
  OutputStream.precision(8);
  
  if(!OutputStream){
        cerr << output_filename << "could not be opened" << endl;
        return;
  }

  for(int i = 0; i < m_grid; i++){
    OutputStream << m_intensity[i] << endl;
  }
  OutputStream.close();
}

    

#endif
