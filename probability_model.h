#ifndef PROBABILITY_MODEL_H
#define PROBABILITY_MODEL_H

#include <gsl/gsl_sf_gamma.h>
#include "gsl/gsl_cdf.h"
#include "changepoint.h"
#include "Data.h"
#include "step_function.h"
#include <algorithm>
#include "particle.h"
#define LOG_TWO log(2.0)

class probability_model{

  public:

  probability_model(Data<double> * m_data, Step_Function* seasonal_scale = NULL);
  probability_model(vector<string>* data_filenames, double start = 0, double end = DBL_MAX, double season = DBL_MAX, bool make_time_scale = true );
  probability_model(string* data_filename = NULL, string* seasonal_data_filename = NULL, double season = DBL_MAX, bool make_time_scale = true );
  virtual ~probability_model();
  void construct();
  virtual double log_likelihood_interval(changepoint *, changepoint *) = 0 ;
  virtual double log_likelihood_interval(double t1, double t2){ return 0;}
  virtual void propose_new_parameters(Particle<changepoint>*, int, unsigned int,changepoint *, changepoint *){};//if third argument 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double calculate_prior_ratio(Particle<changepoint>*,unsigned int){return 0;};//if argument 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double proposal_ratio(Particle<changepoint>*,unsigned int){return 0;};//if 0 birth if 1 death if 2 move changepoint if 3 move parameter
  virtual double calculate_log_predictive_df(double t1, double t2, double t3, bool lower_tail = true ){ return 0; }
  virtual void calculate_sequential_log_predictive_dfs(double start, double end, double increment, bool lower_tail = true, bool two_sided = false, double control_chart_weight = 0.05, string* filename_ptr = NULL, vector<double>* dfs = NULL ){if(dfs) for(unsigned int i=0; i<(end-start)/increment; i++) dfs->push_back(-LOG_TWO);}
  virtual double calculate_log_predictive_df_bounds( double dt, bool lower_tail = true, bool two_sided = false, bool increment_parameters = true ){ m_pvalue_pair = make_pair(-LOG_TWO,-LOG_TWO); return -LOG_TWO; }
  virtual void calculate_posterior_mean_parameters(changepoint *, changepoint *){}
  virtual double calculate_mean(changepoint *, changepoint *) = 0;
  virtual double draw_mean_from_posterior(changepoint *, changepoint *){ return m_mean;}
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
};

bool probability_model::m_seasonal_analysis = false;
bool probability_model::m_p_value_alternative_style = false;
bool probability_model::m_monte_carlo_p_values = false;
unsigned int probability_model::m_p_value_monte_carlo_samples = 100000;
double probability_model::m_current_t = 0;

probability_model::probability_model(vector<string>* data_filenames, double start, double end, double season, bool make_time_scale )
:m_start(start),m_end(end),m_season(season)
{
   construct();
   m_owner_of_data = true;
   m_owner_of_time_scale = m_owner_of_seasonal_scale = false;
   m_data_cont = (data_filenames && data_filenames->size()>0) ? new Data<double>((*data_filenames)[0].c_str()) : NULL;

   if(data_filenames && data_filenames->size()>2){
     m_owner_of_seasonal_scale = true;
     m_seasonal_scale = new Step_Function((*data_filenames)[2],season);
     if( !m_seasonal_scale->get_num_changepoints()){
       delete m_seasonal_scale;
       m_seasonal_scale = NULL;
     }
   }
   if(make_time_scale)
     construct_time_scale(data_filenames,season);
}

probability_model::probability_model(string* data_filename, string* seasonal_data_filename, double season, bool make_time_scale){
   construct();
   m_owner_of_data = true;
   m_owner_of_time_scale = m_owner_of_seasonal_scale = false;
   m_data_cont = data_filename ? new Data<double>(*data_filename) : NULL;
   if(seasonal_data_filename){
     m_owner_of_seasonal_scale = true;
     m_seasonal_scale = new Step_Function(*seasonal_data_filename,season);
     if( !m_seasonal_scale->get_num_changepoints()){
       delete m_seasonal_scale;
       m_seasonal_scale = NULL;
     }
   }
}

probability_model::probability_model(Data<double> * data, Step_Function* seasonal_scale)
:m_data_cont(data)
{
   m_owner_of_data = false;
   m_owner_of_time_scale = m_owner_of_seasonal_scale = false;
   construct();
   m_seasonal_scale = seasonal_scale;
}

probability_model::~probability_model(){
  if(m_data_cont && m_owner_of_data)
    delete m_data_cont;
  bool time_and_season_scale_equal = m_seasonal_scale == m_time_scale;
  if(m_seasonal_scale && m_owner_of_seasonal_scale)
    delete m_seasonal_scale;
  if(m_time_scale && m_owner_of_time_scale && !time_and_season_scale_equal)
    delete m_time_scale;
  if(m_data_seasons)
    delete [] m_data_seasons;
  if(m_rng)
    gsl_rng_free(m_rng);
}

void probability_model::construct(){
  m_time_scale = NULL;
  m_seasonal_scale = NULL;
  m_data_seasons = NULL;
  m_num_data_seasons = 0;
  m_currently_observable = true;
  m_rng = NULL;
  m_pvalue_pair_on_log_scale = false;
}

void probability_model::construct_time_scale(vector<string>* data_filenames, double season){
  m_time_scale = NULL;
  if(data_filenames && data_filenames->size()>1){
    m_time_scale = new Step_Function((*data_filenames)[1],m_end);
    if( !(m_time_scale->get_num_changepoints())){
      delete m_time_scale;
      m_time_scale = m_seasonal_scale;//NULL;
    }
    else if(m_seasonal_scale){
      Step_Function* temp = m_time_scale;
      m_time_scale = temp->multiply_step_function( m_seasonal_scale );
      delete temp;
      m_owner_of_time_scale = true;
    }
    if(m_time_scale && data_filenames->size()>3){
      Step_Function* time_scale2 = new Step_Function((*data_filenames)[3]);
      if( time_scale2->get_num_changepoints()){
	Step_Function* temp = m_time_scale;
	m_time_scale = temp->multiply_step_function( time_scale2 );
	if(temp != m_seasonal_scale)
	  delete temp;
      }
      delete time_scale2;
    }
  }
}

void probability_model::set_data_index(changepoint * cpobj, unsigned int i, changepoint * cpobj_left, changepoint * cpobj_right)
{
  double theta = cpobj->getchangepoint();
  if(m_data_cont){
    unsigned long long int dataindex_left, dataindex_right;
    if(cpobj_left)
      dataindex_left = cpobj_left->getdataindex();
    else dataindex_left = 0;
    if(cpobj_right)
      dataindex_right = cpobj_right->getdataindex()-1;
    else
      dataindex_right = m_data_cont->get_cols()-1;
    cpobj->setdataindex(m_data_cont->find_data_index(theta,i,dataindex_left,dataindex_right));
  }else{//it's a discrete model
    cpobj->setdataindex(static_cast<unsigned long long int>(ceil(theta)));
  }
}

void probability_model::find_data_seasons(){
  m_num_data_seasons = m_seasonal_scale->get_num_changepoints();
  double* cps = m_seasonal_scale->get_changepoints();
  unsigned int num_cols = m_data_cont ? m_data_cont->get_cols() : 0;
  m_data_seasons = new unsigned int[ num_cols ];
  for(unsigned int i = 0; i < num_cols; i++){
    double t = (*m_data_cont)[0][i];
    t -= m_season * static_cast<unsigned int>(t/m_season);
    unsigned int j = 0;
    while(j<m_num_data_seasons && t>cps[j] )
      j++;
    m_data_seasons[i] = j-1;
  }
}

void probability_model::collapse_to_seasons(){
  if(m_data_cont){
    m_data_cont->replace_with_modulo(m_season);
    m_data_cont->sort();
  }
  m_start = 0;
  m_end = m_season;
  collapse_to_seasons_implementation();
}

double probability_model::get_two_sided_discrete_p_value( unsigned int x ){
  double g_x = m_df_values[ x ];
  double df_sum = 0, df_sum2 = 0;
  unsigned int i = 0;
  while( g_x<m_df_values[ i ] || i<=x ){
    if(g_x<m_df_values[ i ])
      df_sum += m_pdf_values[ i ];
    if(i==x || g_x<=m_df_values[ i ])
      df_sum2 += m_pdf_values[ i ];
    i++;
  }
  m_predictive_two_sided_df = 1-df_sum;
  m_predictive_two_sided_df2 = 1-df_sum2;
  m_log_predictive_df = log(m_predictive_two_sided_df+m_predictive_two_sided_df2) - LOG_TWO;
  return m_log_predictive_df;
  if(m_log_predictive_df>-LOG_TWO)
//      m_log_predictive_df = log(1-exp(m_log_predictive_df));
    m_log_predictive_df = m_log_predictive_df + log(exp(-m_log_predictive_df)-1);
  if(m_log_predictive_df<-DBL_MAX)
    m_log_predictive_df =  m_log_pdf_const;
  m_log_predictive_df += LOG_TWO;
  m_log_predictive_df -= log(1+(2*m_survivor_midpoint-1)*(2*m_survivor_midpoint-1));
  return m_log_predictive_df;
}

double probability_model::combine_p_values_from_endpoints(){
  return combine_p_values_from_endpoints( m_monte_carlo_p_values, m_p_value_monte_carlo_samples );
}

double probability_model::combine_p_values_from_endpoints( bool monte_carlo, unsigned int monte_carlo_size ){
  unsigned int num_pvals = m_p_value_endpoints.size();
  if(!num_pvals)
    return -LOG_TWO;
  if(num_pvals==1){
    if(m_p_value_endpoints_log_scale[0])
      return m_p_value_endpoints[0].second + log(1+exp(m_p_value_endpoints[0].first-m_p_value_endpoints[0].second)) - LOG_TWO;
    return log(m_p_value_endpoints[0].first+m_p_value_endpoints[0].second) - LOG_TWO;
    //    double a = log_scale ? exp(iter->first) : iter->first;
    //    double b = log_scale ? exp(iter->second) : iter->second;
    //    return log(a+b) - LOG_TWO;
  }
  double combined_value = 0;
  if( monte_carlo && !m_rng){
    m_rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set (m_rng,0);
  }
  unsigned int counter = 0;
  double sum_widths = 0;
  do{
    double sum_log_pvals = 0;
    for(unsigned int i = 0; i < num_pvals; i++){
      if(!monte_carlo){
	if(m_p_value_endpoints_log_scale[i])
	  sum_log_pvals += m_p_value_endpoints[i].second + log(1+exp(m_p_value_endpoints[i].first-m_p_value_endpoints[i].second)) - LOG_TWO;
	else
	  sum_log_pvals += log(m_p_value_endpoints[i].first+m_p_value_endpoints[i].second) - LOG_TWO;
	//	sum_log_pvals += log(a+b) - LOG_TWO;
      }else{
	double a = m_p_value_endpoints_log_scale[i] ? exp(m_p_value_endpoints[i].first) : m_p_value_endpoints[i].first;
	double b = m_p_value_endpoints_log_scale[i]  ? exp(m_p_value_endpoints[i].second) : m_p_value_endpoints[i].second;
	sum_log_pvals += log(gsl_ran_flat(m_rng,a,b));
	if(!counter)
	  sum_widths += fabs(b-a);
      }
    }
    double fisher_score = gsl_cdf_chisq_Q(-2*sum_log_pvals,2*num_pvals);
    combined_value += fisher_score;
    //    double stouffer_score = gsl_cdf_ugaussian_Qinv(fisher_score);
    //    combined_value += stouffer_score;
    if(monte_carlo && sum_widths<=0)
      monte_carlo_size = 1;
  } while(monte_carlo && ++counter < monte_carlo_size );
  if(monte_carlo)
    combined_value /= monte_carlo_size;
  //  combined_value = gsl_cdf_ugaussian_Q(combined_value);
  return log(combined_value);

  /*  if(!monte_carlo){
    for(vector< pair<double,double> >::iterator iter = m_p_value_endpoints.begin(); iter != m_p_value_endpoints.end(); ++iter){
      if(log_scale)
	sum_log_pvals += log(1+exp(iter->second-iter->first)) + iter->first - LOG_TWO;
      else
	sum_log_pvals += log(iter->second+iter->first) - LOG_TWO;
    }
    combined_value = gsl_cdf_chisq_Q(-2*sum_log_pvals,2*m_p_value_endpoints.size());
  }else{
    for(unsigned int i = 0; i<monte_carlo_size; i++){
      sum_log_pvals = 0;
      for(vector< pair<double,double> >::iterator iter = m_p_value_endpoints.begin(); iter != m_p_value_endpoints.end(); ++iter)
	sum_log_pvals += log(gsl_ran_flat(m_rng,log_scale?exp(iter->first):iter->first,log_scale?exp(iter->second):iter->second));
      combined_value += gsl_cdf_chisq_Q(-2*sum_log_pvals,2*m_p_value_endpoints.size());
    }
    combined_value/=monte_carlo_size;
  }
  return log(combined_value);*/
}


#endif
