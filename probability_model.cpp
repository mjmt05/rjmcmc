#include "probability_model.hpp"

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
  m_num_windows = 0;
  m_windows = NULL;
  m_windowed_lhd_contributions = NULL;
  window_mixture_probs = NULL;
  window_mixture_probs_sum = 0;
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
