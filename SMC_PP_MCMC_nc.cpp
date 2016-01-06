#include "SMC_PP_MCMC_nc.hpp"
#define LOG_TWO log(2.0)
#include "string.h"
#include <iostream>

bool SMC_PP_MCMC::MyDataSort(const pair<double,int>& lhs, const pair<double,int>& rhs){return (lhs.first > rhs.first);}

SMC_PP_MCMC::SMC_PP_MCMC(double start, double end, unsigned int intervals, int sizeA, int sizeB, double nu, double v_nu, probability_model ** pm,int num_data,bool varyB,bool intensity,bool dochangepoint,bool doMCMC, bool exact_sampling, int s)
  :SMC_PP<changepoint>(start,end,intervals,sizeA,sizeB,num_data,varyB,dochangepoint,doMCMC,s),m_calculate_intensity(intensity), m_do_exact_sampling(exact_sampling)
{

  m_pm = pm;
  m_nu = nu;
  m_var_nu = v_nu;

  m_prior_diff=0;
  m_proposal_prior=m_nu;

  m_conjugate=true;
  m_use_spacing_prior = false;
  m_functionofinterest=NULL;
  m_do_SMC_past=0;
  m_length_grid=0;
  if (m_do_exact_sampling) {
    m_variable_B = 0;
  }
  m_rejection_sampling_acceptance_rate = NULL;
  m_num_zero_weights = NULL;

  if(MCMC_only){
    m_rj_A=new rj_pp*[m_num];
    for(int ds=0; ds<m_num; ds++){
      m_rj_A[ds]=NULL;}
    m_rj_B=NULL;
    m_rejection_sampling = NULL;
    m_do_exact_sampling = 0;
  }else{
    if (m_do_exact_sampling) {
      m_rejection_sampling = new rejection_sampling *[m_num];
      for (int ds = 0; ds < m_num; ds++) {
	m_rejection_sampling[ds] = NULL;
      }
      m_rejection_sampling_acceptance_rate = new double *[m_num];
      m_rejection_sampling_acceptance_rate[0] = new double[m_num * m_num_of_intervals];
      for (int i = 1; i < m_num; i++) {
	m_rejection_sampling_acceptance_rate[i] = m_rejection_sampling_acceptance_rate[i-1] + m_num_of_intervals;
      }
    } else {
      m_rj_B=new rj_pp*[m_num];
      for(int ds=0; ds<m_num; ds++){
	m_rj_B[ds]=NULL;
      }
    }
    m_rj_A=NULL;
  }

  m_spacing_prior = m_change_in_time;
  
  m_neighbouring_intervals=0;

  m_calculate_intensity=intensity;
 
  m_move_width=(m_end-m_start)/10.0;
  m_thin=10;
  if(!m_variable_B)
    m_burnin=(int)(m_max_sample_size_A/100);
  else
    m_burnin=1000;

  m_store_sample_sizes=0;
  m_sample_sizes=NULL;
  if(m_variable_B){
    m_min_sample_size = new unsigned long long int[m_num]; 
    m_current_sample_size = new unsigned long long int[m_num];
    for(int i=0; i<m_num; i++) {
      m_min_sample_size[i]=m_max_sample_size_A;
      m_current_sample_size[i] = m_max_sample_size_A / m_num;
    }
    m_vec_KLS = new double[m_num];
   
  }else{
    m_foi_grid=0;}

  m_proposal_type=NULL;

}


SMC_PP_MCMC::~SMC_PP_MCMC(){
  
  for(int ds=0; ds<m_num; ds++){
    if(MCMC_only){
      delete m_rj_A[ds];
    }else{
      if (m_rj_B) {
	delete m_rj_B[ds];
      } else {
	delete m_rejection_sampling[ds];
      }
      delete_samples(ds);
      delete [] m_sample_A[ds];
    }
  }

  if(m_functionofinterest){
    for(int ds=0; ds<m_num; ds++)     
      delete m_functionofinterest[ds];
    delete [] m_functionofinterest;
  }

 
  if(MCMC_only){
    delete [] m_rj_A;}
  else{
    if (m_rj_B) {
      delete [] m_rj_B;
    } else {
      delete [] m_rejection_sampling_acceptance_rate[0];
      delete [] m_rejection_sampling_acceptance_rate;
      delete [] m_rejection_sampling;
    }
  }
  
  if (m_num_zero_weights) {
    delete [] m_num_zero_weights[0];
    delete [] m_num_zero_weights;
  }

  if(m_variable_B){
    delete [] m_current_sample_size;
    delete [] m_vec_KLS;
    delete [] m_min_sample_size;
    if(m_store_sample_sizes){
      delete [] m_sample_sizes[0];
      delete [] m_sample_sizes;
    }
  }

}

void SMC_PP_MCMC::set_RJ_parameters(int t, int b, double mw, const char *proposaltype, void * v){
  m_thin=t;
  m_burnin=b; 
  m_move_width=mw;

  if(proposaltype){
    const char *no = "Normal";
    const char *mix = "Mixture";
    const char *hist = "Histogram";

    m_vec_proposal_type=v;
    if(strcmp(no,proposaltype) == 0){
      m_proposal_type = no;

    }
    else if(strcmp(mix,proposaltype)==0){
      m_proposal_type=mix;
  
    }else if (strcmp(hist,proposaltype)==0){
      m_proposal_type = hist;
    }
  }
}

void SMC_PP_MCMC::use_spacing_prior() {
  m_use_spacing_prior = true;
  m_num_zero_weights = new unsigned int *[m_num];
  m_num_zero_weights[0] = new unsigned int[m_num * m_num_of_intervals];
  for (int i = 1; i < m_num; i++) {
    m_num_zero_weights[i] = m_num_zero_weights[i-1] + m_num_of_intervals;
  }
  for (int ds = 0; ds < m_num; ds++) {
    for (unsigned int i = 0; i < m_num_of_intervals; i++) {
      m_num_zero_weights[ds][i] = 0;
    }
  }
}

void SMC_PP_MCMC::initialise_function_of_interest(int grid, bool g, bool prob, bool instant, double delta, bool sequential){

  m_length_grid=grid;
  m_functionofinterest = new Function_of_Interest * [m_num];
  for(int ds=0; ds<m_num; ds++){
    m_functionofinterest[ds] = new Function_of_Interest(m_length_grid,m_start,m_end,m_nu,g,prob,m_calculate_intensity,1,sequential,m_num_of_intervals,instant,delta);
    
    if (m_importance_sampling) {
      m_functionofinterest[ds]->set_importance_sampling();
    }
  }
    
}

void SMC_PP_MCMC::delete_samples(int ds){
    
  for(unsigned long long int i=m_sample_size_A[ds]-1; i>0; i--){
    if(m_sample_A[ds][i]!=m_sample_A[ds][i-1]){
      delete m_sample_A[ds][i];
      
    }
  }
  delete m_sample_A[ds][0];

  if(m_sample_size_B[ds]>m_sample_size_A[ds]){
    m_sample_size_A[ds]=m_sample_size_B[ds];
  }
   
}

void SMC_PP_MCMC::store_sample_sizes(){
  m_store_sample_sizes=1;
  m_sample_sizes = new unsigned long long int * [m_num];
  m_sample_sizes[0] = new unsigned long long int[m_num*m_num_of_intervals];
  for(int i=1; i<m_num; i++)
    m_sample_sizes[i]=m_sample_sizes[i-1]+m_num_of_intervals;

  for(int i=0; i<m_num; i++){
    for(unsigned int j=0; j<m_num_of_intervals;j++)
      m_sample_sizes[i][j]=0;
  }
}

void SMC_PP_MCMC::sample_particles(double start, double end){
  
  long double avg_distance=0;
  //  double variance=0;
  double move_width=m_move_width;
  bool active=0;
  unsigned long long int current_number=0; 
  
  Particle<changepoint> * tempparticle;      
  for(int ds=0; ds<m_num; ds++){
    if(MCMC_only==1){
      if(iters>0){
	tempparticle = new Particle<changepoint>(m_rj_A[ds]->get_MAP_dimension_MAP(),NULL);
	changepoint * cpobj1 = new changepoint(end,0,0,0);
	m_pm[ds]->set_data_index(cpobj1);
	int position=tempparticle->get_dim_theta()-1;
	if(!m_conjugate)
	  m_pm[ds]->propose_new_parameters(tempparticle,position,3,NULL,cpobj1);
	changepoint * cpobj = tempparticle->get_theta_component(position);
	cpobj->setlikelihood(m_pm[ds]->log_likelihood_interval(cpobj,cpobj1,position>=0?tempparticle->get_theta_component(position):NULL));
	cpobj->setmeanvalue(m_pm[ds]->calculate_mean(cpobj, cpobj1,position>=0?tempparticle->get_theta_component(position):NULL));
	delete m_rj_A[ds];
	delete cpobj1;
      }
      else{
	tempparticle=NULL;
      }
      m_rj_A[ds] = new rj_pp(start, end, m_sample_size_A[ds], 100000, move_width , m_nu, m_var_nu, m_pm[ds],m_thin,m_burnin,0,0,tempparticle,(int)(seed*(iters+1)),true);
      if(m_proposal_type && m_vec_proposal_type){
	if(strcmp(m_proposal_type,"Histogram")==0){
	  int temp=(((unsigned int*)m_vec_proposal_type)[0])*(iters+1);
	  m_rj_A[ds]->proposal_type(m_proposal_type,(void*)(&temp));
	}
      }
     
      if(!m_conjugate)
	m_rj_A[ds]->non_conjugate();
      if(m_calculate_intensity && !m_importance_sampling)
	m_rj_A[ds]->calculate_intensity();
      if(m_neighbouring_intervals)
	m_rj_A[ds]->allow_neighbouring_empty_intervals();
      else
	m_rj_A[ds]->disallow_neighbouring_empty_intervals();
     
      m_rj_A[ds]->runsimulation();
      m_sample_A[ds] = m_rj_A[ds]->get_sample();
      if (m_importance_sampling) {
	sample_intensities(m_sample_A[ds], end, m_sample_size_A[ds], ds);
      }
    }else{  

      if(m_process_observed[ds]==0){
	active=1;
	m_process_observed[ds]++;
      }else{      
	active=1;
	if (!m_do_exact_sampling) {
	  delete m_rj_B[ds];
	} else {
	  delete m_rejection_sampling[ds];
	}
	m_process_observed[ds]++;
      }

      if(m_process_observed[ds]>0){
	if(m_do_SMC_past && m_process_observed[ds]>1){
	  avg_distance=m_functionofinterest[ds]->get_average_distance();
	  //  variance=m_functionofinterest[ds]->get_variance_distance();
	}else if(m_process_observed[ds]==1){
	    avg_distance=start-m_start;
	}else{
	  avg_distance = 0;
	}

	unsigned long long int sample_size=0;

	sample_size = m_max_sample_size_A;
	if(!m_variable_B && m_process_observed[ds]>1){
	  sample_size=m_max_sample_size_B;
	}
	
	unsigned int max_theta = UINT_MAX;
	m_cp_start = start;
	if (m_use_spacing_prior && !m_sample_from_prior) {
	  max_theta = 1;
	  if (start > m_start) {
	    m_cp_start = max(m_functionofinterest[ds]->get_min_distance() + m_spacing_prior, start);
	  }
	}
	if (!m_do_exact_sampling) {
	  m_rj_B[ds] = new rj_pp((double)(start-avg_distance), (double)end, sample_size, max_theta ,move_width, m_nu, m_var_nu, m_pm[ds],m_thin,m_burnin,0,m_variable_B,NULL,(int)seed*(iters+1),true);
	  
	  if(!m_conjugate){
	    m_rj_B[ds]->non_conjugate();
	  }
	  if(m_calculate_intensity && !m_importance_sampling){
	    m_rj_B[ds]->calculate_intensity();
	  }
	  m_rj_B[ds]->set_start_cps(m_cp_start);
	  if (m_use_spacing_prior) {
	    m_rj_B[ds]->set_spacing_prior();
	  }
	  if(m_proposal_type && m_vec_proposal_type){
	    m_rj_B[ds]->proposal_type(m_proposal_type,m_vec_proposal_type);
	  }
	  if(m_neighbouring_intervals){
	    m_rj_B[ds]->allow_neighbouring_empty_intervals();
	  }else{
	    m_rj_B[ds]->disallow_neighbouring_empty_intervals();
	  }

	  if(m_variable_B){
	    m_rj_B[ds]->set_initial_iterations(m_initial_iterations);
	    if(m_divergence_type==FOI){
	      m_rj_B[ds]->set_function_criteria(m_foi_grid);
	    }
	    m_rj_B[ds]->start_calculating_mc_divergence(m_divergence_type,m_loss_type,m_num_of_bins,true,(unsigned int)(m_max_sample_size_A<m_max_lookup_length?m_max_sample_size_A:m_max_lookup_length));
	    m_rj_B[ds]->no_1d_histogram();
	  }

	  m_rj_B[ds]->runsimulation();

	  if(m_variable_B){
	    m_vec_KLS[ds]=m_rj_B[ds]->get_divergence();
	    current_number+=m_initial_iterations;
	    m_rj_B[ds]->end_divergence_burn_in();
	    m_rj_B[ds]->set_continue_loop(1);
	  }    
	  m_sample_B[ds] = m_rj_B[ds]->get_sample();
	  if (m_importance_sampling) {
	    sample_intensities(m_sample_B[ds], end, sample_size, ds);
	  }
	} else {
	 
	 
	  m_rejection_sampling[ds] = new rejection_sampling((double)(start - avg_distance),  m_cp_start, (double) end, 
							    sample_size, m_pm[ds], m_nu, (int) seed * (iters + 1),
							    m_calculate_intensity);
	  if (!m_sample_from_prior) {
	    m_rejection_sampling[ds]->run_simulation();
	    m_rejection_sampling_acceptance_rate[ds][iters] = m_rejection_sampling[ds]->m_acceptance_rate;
	  } else {
	    if (start == m_start) {
	      m_rejection_sampling[ds]->sample_from_prior(NULL);
	    } else {
	      m_rejection_sampling[ds]->sample_from_prior(m_sample_A[ds]);
	    }
	  }
	  m_sample_B[ds] = m_rejection_sampling[ds]->get_sample();
	  
	  //cerr << m_rejection_sampling[ds]->m_acceptance_rate << endl;
	 
	  
	}	  
      } else {
	if(m_variable_B){
	  m_vec_KLS[ds]=0;
	}
      }
    }
  }

  if(m_variable_B && active && !MCMC_only){
    int current_max = find_max(m_vec_KLS,m_num);
   
    while(current_number<m_max_sample_size_A){
         m_rj_B[current_max]->runsimulation();
      m_vec_KLS[current_max]=m_rj_B[current_max]->get_divergence();
      m_rj_B[current_max]->set_continue_loop(1);
      current_number+=1;
      current_max=find_max(m_vec_KLS,m_num);
    }
    for(int ds=0; ds<m_num; ds++){
      if(m_process_observed[ds]==1){
	if (m_rj_B[ds]->get_size_sample() > m_current_sample_size[ds]) {
	  increase_vector(ds);
	}  
	m_sample_size_A[ds]=m_rj_B[ds]->get_size_sample();
	if(m_store_sample_sizes)
	  m_sample_sizes[ds][iters]=m_sample_size_A[ds];
	if(m_sample_size_A[ds]<m_min_sample_size[ds])
	  m_min_sample_size[ds]=m_sample_size_A[ds];
      }else if(m_process_observed[ds]>1){
	if (m_rj_B[ds]->get_size_sample() > m_current_sample_size[ds]) {
	  increase_vector(ds);
	}
	m_sample_size_B[ds]=m_rj_B[ds]->get_size_sample();
	if(m_store_sample_sizes)
	  m_sample_sizes[ds][iters]=m_sample_size_B[ds];
	if(m_sample_size_B[ds]<m_min_sample_size[ds])
	  m_min_sample_size[ds]=m_sample_size_B[ds];
      }else
	{m_sample_size_A[ds]=0; m_sample_size_B[ds]=0;}
    }
  }
}

void SMC_PP_MCMC::sample_intensities(Particle<changepoint> ** sample, double end, unsigned int sample_size, int ds) {
  int dim;
  changepoint *cp;
  changepoint *cp1;
  changepoint *cp0;
  changepoint *end_of_int_changepoint = new changepoint(end, 0, 0, 0);
  m_pm[ds]->set_data_index(end_of_int_changepoint);
  for (unsigned int i = 0; i < sample_size; i++) {
    dim = sample[i]->get_dim_theta();
    cp = sample[i]->get_theta_component(-1);
    cp0 = NULL;
    for (int j = 0; j < dim; j++) {
      cp1 = sample[i]->get_theta_component(j);
      cp->setmeanvalue(m_pm[ds]->calculate_mean(cp, cp1, cp0));
      cp0 = cp;
      cp = cp1;
    }
    cp->setmeanvalue(m_pm[ds]->calculate_mean(cp, end_of_int_changepoint,cp0));
  }
}

void SMC_PP_MCMC::increase_vector(int ds) {
  unsigned long long int size = m_current_sample_size[ds] * 2;
  unsigned long long int j;
  m_current_sample_size[ds] = size;
  if (m_process_observed[ds] > 1) {
    for (j = 0; j < m_sample_size_A[ds]; j++) {
      m_sample_dummy[ds][j] = m_sample_A[ds][j];
    }
  }
  delete [] m_sample_A[ds];
  m_sample_A[ds] = new Particle <changepoint> *[size];
  
  if (m_process_observed[ds] > 1) {
    for (j = 0; j < m_sample_size_A[ds]; j++) {
      m_sample_A[ds][j] = m_sample_dummy[ds][j];
    }
  }
   
  delete [] m_sample_dummy[ds];
  m_sample_dummy[ds] = new Particle <changepoint> *[size];
  delete [] m_exp_weights[ds];
  delete [] m_cum_exp_weights[ds];
  m_exp_weights[ds] = new double[size];
  m_cum_exp_weights[ds] = new double[size];
  for (j = 0; j < size; j++) {
    m_exp_weights[ds][j] = 1;
  }
  double *weightstemp = NULL;
  if (m_process_observed[ds] > 1) {
    weightstemp = new double[m_sample_size_A[ds]];
    for (unsigned int i = 0; i < m_sample_size_A[ds]; i++) {
      weightstemp[i] = m_weights[ds][i];
    }
  }
  delete [] m_weights[ds];
  m_weights[ds] = new double[size];

  if (m_process_observed[ds] > 1) {
    for (unsigned int i = 0; i < m_sample_size_A[ds]; i++) {
      m_weights[ds][i] = weightstemp[i];
    }
    delete [] weightstemp;
  }
}

void SMC_PP_MCMC::resample_particles(double start, double end, int num, const char * ptr2char,int ds){
    double move_width=m_move_width;
    double normal_pars[2] = {start,(end-start)/3};

    rj_pp * rj_pp_obj = new rj_pp(start,end,num,10000,move_width,m_nu,m_var_nu,m_pm[ds],1,0,0,0,NULL,seed*(iters+1),true);

    if(m_proposal_type && m_vec_proposal_type){
      if(strcmp(m_proposal_type,"Histogram")==0){
	int temp=(((unsigned int*)m_vec_proposal_type)[0])*(iters+1);
	rj_pp_obj->proposal_type(m_proposal_type,(void*)(&temp));
      }else
	rj_pp_obj->proposal_type(m_proposal_type,m_vec_proposal_type);
    }else if(strcmp(ptr2char,"Uniform")!=0){
      rj_pp_obj->proposal_type(ptr2char,(void*)(&normal_pars));
    }

    if(m_calculate_intensity){
      rj_pp_obj->calculate_intensity();
    }

    if(!m_conjugate){
      rj_pp_obj->non_conjugate();
    }

    rj_pp_obj->stop_storing_sample();

    for(unsigned long long int i=m_sample_size_A[ds]-1; i>0; i--){
        if(m_sample_A[ds][i]==m_sample_A[ds][i-1]){
            m_sample_A[ds][i]=new Particle<changepoint>(m_sample_A[ds][i],NULL);
	}
        rj_pp_obj->set_initial_sample(m_sample_A[ds][i]);
       	rj_pp_obj->runsimulation();
    }
    
    rj_pp_obj->set_initial_sample(m_sample_A[ds][0]);
    rj_pp_obj->runsimulation();
    m_sample_A[ds][0] = new Particle<changepoint>(rj_pp_obj->get_current_particle(),NULL);
    delete rj_pp_obj;
}

void SMC_PP_MCMC::ESS_resample_particles(double end,int ds){

    int * num_resampled_particles = new int[ m_sample_size_A[ds]];
    double unif_rand = gsl_ran_flat(r,0,1)*(1.0/((double)m_sample_size_A[ds]));
 
    for (unsigned int i=0; i<m_sample_size_A[ds]; i++){
        num_resampled_particles[i]=0;
    }

    unsigned int  j=0, k=0;

    while(k<m_sample_size_A[ds]){
      while((m_cum_exp_weights[ds][j]/m_sum_exp_weights[ds] - unif_rand) > ((double)k)/((double)m_sample_size_A[ds]) && k<m_sample_size_A[ds]){
	num_resampled_particles[j]++;
	k++;
      }
      j++;
    }

    Particle<changepoint> ** holding;
    holding = new Particle<changepoint> *[m_max_sample_size_A];

    for(unsigned int i=0; i<m_sample_size_A[ds]; i++){
        if(num_resampled_particles[i]==0)
            delete m_sample_A[ds][i];
    }
        
    j=0;

    for(unsigned int i=0; i<m_sample_size_A[ds]; i++){
      while(num_resampled_particles[i]>0){
	holding[j++]=m_sample_A[ds][i];
	-- num_resampled_particles[i];
      }
    }

    delete [] m_sample_A[ds];
    m_sample_A[ds] = holding;
    delete [] num_resampled_particles;

    for(unsigned int i=0; i<m_sample_size_A[ds]; i++){
      m_weights[ds][i]=0;
    }
}

void SMC_PP_MCMC::increase_A_particles(int ds, unsigned long long int size_increase,unsigned long long int * m){

  list<pair<double, int> > delta;  
  list<pair<double, int> >::iterator iter;
  pair<double,int> p;
  double term;
  unsigned long long int n0=0;
  vector<int> zero_index;

  for(unsigned int i=0; i<m_sample_size_A[ds]; i++){
    if(m_sample_A[ds][i]->get_dim_theta()>0){
      term = 2*m_weights[ds][i]-LOG_TWO;//-log(1);
      delta.push_back(make_pair(term,i));
    }else{
      zero_index.push_back(i);}
  }

  n0=zero_index.size();
  if(n0<m_sample_size_A[ds]){

    for(unsigned int i=0; i<n0; i++){
      m_weights[ds][zero_index[i]]+=log(n0);
      if(i==0){
	term=2*(m_weights[ds][zero_index[i]])-log(n0)-log(n0+1);
	delta.push_back(make_pair(term,zero_index[i]));
	m[zero_index[i]]=n0;
      }
    }
    

    delta.sort(MyDataSort);
    
    unsigned long long int how_many=0;
    unsigned long long int number_increase = 0;
    int index =0;
    bool c;
 
    while(number_increase<size_increase){

      if(number_increase>0){
	p = make_pair(2*m_weights[ds][index]-log(m[index])-log(m[index]+1),index);
	delta.pop_front();
	c=0;
	iter=delta.begin();
	while(!c && iter!=delta.end()){
	  c=MyDataSort(p,*iter);
	  if(!c)
	    ++iter;
	}
	if(iter==delta.begin()){
	  cerr<<"SMC_PP_MCMC_nc: problem"<<endl;
	  exit(1);
	}
	delta.insert(iter,p);

      }
      iter=delta.begin();
      index = (*iter).second;
      term=(*iter).first;
      ++iter;

      if((*iter).first==term){
	how_many=1;
      }else{
	how_many=(unsigned long long int)ceil(-0.5-(double)m[index]+sqrt(exp(2*m_weights[ds][index]-(*iter).first)+0.25));}
       
      if(how_many==0){
	how_many=1;
      }
    
      if(how_many>(size_increase-number_increase)){
    	how_many=size_increase-number_increase;
      }
  
      m[index]+=how_many;
      number_increase+=how_many;

    }

    for(unsigned int i=0; i<(m_sample_size_A[ds]); i++){
      if(m_sample_A[ds][i]->get_dim_theta()>0){
	m_weights[ds][i]-=log(m[i]);
      }
    }

    if(n0>0){
      int zi = zero_index[0];
      for(unsigned int i=0; i<zero_index.size(); i++){
	m_weights[ds][zero_index[i]]-=log(m[zi]); 
      }
      m[zi]-=n0-1;
    }
  }else{
    m[0]+=size_increase;
  }
}


void SMC_PP_MCMC::calculate_weights_join_particles(int iter,int ds){

  if (m_process_observed[ds]==1){
    for(unsigned int i=0; i<m_sample_size_A[ds]; i++){
	m_weights[ds][i]=0; 
	m_sum_exp_weights[ds]=(double)(m_sample_size_A[ds]);
	m_sum_squared_exp_weights[ds]=(double)(m_sample_size_A[ds]);
	m_sample_dummy[ds][i] = new Particle<changepoint>(m_sample_B[ds][i],NULL);
	if (m_sample_from_prior) {
	  m_weights[ds][i] += m_sample_dummy[ds][i]->get_theta_component(-1)->getlikelihood();
	  for (unsigned int j = 0; j < m_sample_dummy[ds][i]->get_dim_theta(); j++) {
	    m_weights[ds][i] += m_sample_dummy[ds][i]->get_theta_component(j)->getlikelihood();
	  }
	}
    }
    m_nu = m_proposal_prior;
  } else if(m_process_observed[ds]>1){
    int dim=0;
    int dim1=0;
    long double likelihood_joint=0;long double likelihood_left=0;long double likelihood_right =0;
    changepoint *cpobjA=NULL;
    changepoint *cpobjB=NULL;
    changepoint *cpobjB1=NULL;
    long double incremental_weight;
    unsigned int index;
    unsigned long long int * m_A = new unsigned long long int[m_sample_size_A[ds]];
    unsigned long long int * m_B = new unsigned long long int[m_sample_size_B[ds]];
    long double prior_term=0;
    long double clean_prior_term = 0;
    bool gt = 0;
    bool weights0 = 0;
    double start = 0;
    if (m_use_spacing_prior) {
      start = m_start + m_change_in_time * iter;
      if (m_cp_start > start) {
	gt = 1;
	clean_prior_term = -m_nu * (m_cp_start - m_spacing_prior);
      } else {
	clean_prior_term = m_nu * m_spacing_prior;
      }
    }
    unsigned long long int ratio=m_sample_size_A[ds];
    for(unsigned int i=0; i<m_sample_size_A[ds]; i++)
      m_A[i]=1; 
    for(unsigned int i=0; i<m_sample_size_B[ds]; i++)
      m_B[i]=1;

    if(m_sample_size_A[ds]<m_sample_size_B[ds]){
      ratio = m_sample_size_B[ds]-m_sample_size_A[ds];

      if(iters>1){
	increase_A_particles(ds,ratio,m_A);
      }else{
	index=0;
	while(ratio>0){
	  if(index==(m_sample_size_A[ds])){index=0;}
	  m_A[index]++;
	  index++;
	  ratio--;
	}
      }
      ratio=m_sample_size_B[ds];
    }else if(m_sample_size_A[ds]>m_sample_size_B[ds]){
      ratio=m_sample_size_A[ds]-m_sample_size_B[ds];
      index=0;
      while(ratio>0){
	if(index==(m_sample_size_B[ds])){index=0;}
	m_B[index]++;
	index++;
	ratio--;
      }
      ratio=m_sample_size_A[ds];
    }

    changepoint * cpobjBalt = new changepoint(m_start+m_change_in_time*(iter+1),0,0,0);
    m_pm[ds]->set_data_index(cpobjBalt,0,m_sample_B[ds][0]->get_theta_component(-1));
	
    int index_B=0;     
    int index_A=0;
    unsigned long long int counter_A=m_A[index_A];
    unsigned long long int counter_B=m_B[index_B];
    long double * dummy_weights = new long double[ratio];

    for(unsigned int index_new=0; index_new<ratio; index_new++){
      if(counter_A==m_A[index_A]){
	dim=m_sample_A[ds][index_A]->get_dim_theta();
	cpobjA = m_sample_A[ds][index_A]->get_theta_component(dim-1); 
	likelihood_left = cpobjA->getlikelihood();
	if (m_use_spacing_prior) {
	  if (gt) {
	    prior_term = clean_prior_term - m_nu * cpobjA->getchangepoint();
	  } else {
	    prior_term = clean_prior_term - m_nu * min(m_spacing_prior, start - cpobjA->getchangepoint());
	  }
	}
      }
	
      if(counter_B==m_B[index_B]){
	dim1=m_sample_B[ds][index_B]->get_dim_theta();
	cpobjB = m_sample_B[ds][index_B]->get_theta_component(-1);//intercept
	likelihood_right = cpobjB->getlikelihood();
	if (!m_use_spacing_prior) {
	  prior_term=dim1*m_prior_diff;
	}
	if(dim1==0){
	  cpobjB1 = cpobjBalt;
	}else{                
	  cpobjB1=m_sample_B[ds][index_B]->get_theta_component(0);
	  if (m_use_spacing_prior && 
	      (cpobjA->getchangepoint() + m_spacing_prior > cpobjB1->getchangepoint() 
	       || isinf(m_weights[ds][index_new]))) {
	    weights0 = 1;
	  }
	}
      }

      if (!weights0) {
	incremental_weight = 0;
	if(!m_conjugate){
	  m_pm[ds]->propose_combined_parameters(m_sample_A[ds][index_A],m_sample_B[ds][index_B],cpobjBalt,m_start+m_change_in_time*iter);
	}
	
	m_sample_dummy[ds][index_new] = new Particle<changepoint>(m_sample_A[ds][index_A],m_sample_B[ds][index_B]);
	
	if (!m_sample_from_prior) {
	  changepoint * cpobj_new_A = m_sample_dummy[ds][index_new]->get_theta_component(dim-1);

	  likelihood_joint = m_pm[ds]->log_likelihood_interval(cpobj_new_A,cpobjB1,dim>0?m_sample_dummy[ds][index_new]->get_theta_component(dim-2):NULL);
	
	  cpobj_new_A->setlikelihood(likelihood_joint);
	  double mean;
	  if(m_calculate_intensity && m_conjugate){
	    mean = m_pm[ds]->calculate_mean(cpobj_new_A,cpobjB1,dim>0?m_sample_dummy[ds][index_new]->get_theta_component(dim-2):NULL);
	    cpobj_new_A->setmeanvalue(mean);
	  }
        
	  incremental_weight += likelihood_joint-likelihood_left-likelihood_right+prior_term;
	} else {
	  incremental_weight += m_sample_B[ds][index_B]->get_theta_component(-1)->getlikelihood();
	  for (int i = 0; i < dim1; i++) {
	    incremental_weight += m_sample_B[ds][index_B]->get_theta_component(i)->getlikelihood();
	  }
	}
	  
	double b;
	if(!m_conjugate){
	  //any prior ratios should be calculated in here
	  b=m_pm[ds]->non_conjugate_weight_terms(m_sample_dummy[ds][index_new]);
	  incremental_weight+=b;
	}
	dummy_weights[index_new]=m_weights[ds][index_A]+incremental_weight;
      } else {
	dummy_weights[index_new] = log(0); 
	changepoint *cp = new changepoint(cpobjB);
	m_sample_dummy[ds][index_new] = new Particle<changepoint>(0, NULL, cp);
	incremental_weight = 0;
	weights0 = 0;
	m_num_zero_weights[ds][iter]++;
      }
      
      
      m_sample_dummy[ds][index_new]->set_log_posterior((double)m_sample_dummy[ds][index_new]->get_log_posterior() + (double)incremental_weight);

      if(index_new!=(ratio-1)){
	if(--counter_A==0){	   
	  index_A++;
	  counter_A=m_A[index_A];
	}
	    
	if(--counter_B==0){	   
	  index_B++;
	  counter_B=m_B[index_B];
	}
      }
    }
	
    for(unsigned int i=0; i<ratio; i++){
      //     if(isinf(m_weights[ds][i])) {
   	// }
      m_weights[ds][i]=(double)(dummy_weights[i]);
    }
   
        
    delete cpobjBalt;
    delete [] m_A;
    delete [] m_B;
    delete [] dummy_weights;
  }
}

void SMC_PP_MCMC::calculate_function_of_interest(double start, double end){
  if (m_functionofinterest){
    for(int ds=0; ds<m_num; ds++){
      if(start==m_start && end==m_end){

	m_functionofinterest[ds]->reset_prob();
	}
      if(m_process_observed[ds]>0)
	m_functionofinterest[ds]->calculate_function(start,end,m_sample_A[ds],m_sample_size_A[ds],m_exp_weights[ds],m_sum_exp_weights[ds],m_sum_squared_exp_weights[ds],iters,1,m_pm[ds]);
    }
  }
}

void SMC_PP_MCMC::print_intensity(int ds, const char* file){
   m_functionofinterest[ds]->write_mean_to_file(file);
}


void SMC_PP_MCMC::print_zero_weights(int ds, const char *file) {
  if (!m_num_zero_weights) {
    return;
  }
  
  ofstream outfile(file, ios::out);
  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }

  for (unsigned int i = 0; i < m_num_of_intervals; i++) {
    outfile << m_num_zero_weights[ds][i] << endl;
  }
  outfile.close();
}

void SMC_PP_MCMC::print_rejection_sampling_acceptance_rates(int ds, const char *file) {
  if (!m_rejection_sampling_acceptance_rate) {
    return;
  }
  
  ofstream outfile(file, ios::out);
  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }

  for (unsigned int i = 0; i < m_num_of_intervals; i++) {
    outfile << m_rejection_sampling_acceptance_rate[ds][i] << endl;
  }
  outfile.close();
}


void SMC_PP_MCMC::print_exp(const char* file){
  ofstream outfile(file, ios::out);
  outfile<<setiosflags(ios::fixed);
  outfile.precision(5);

  if(!outfile){
        cerr<<outfile << " could not be opened"<<endl;
	return;
    }

  for(int ds=0; ds<m_num; ds++){
    long  double * g=m_functionofinterest[ds]->get_g();
    
    for(int i=0; i<m_length_grid; i++){
      outfile<<g[i]<<" ";
    }
    outfile<<endl;
  }
  outfile.close();
}


void SMC_PP_MCMC::print_var_exp(int ds,const char* file){

  long  double * g=m_functionofinterest[ds]->get_variance_g();
  ofstream outfile(file, ios::out);

  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }

  for(int i=0; i<m_length_grid; i++){
    outfile<<g[i]<<endl;
  }
  outfile.close();
}

void SMC_PP_MCMC::print_exp_sequential(int ds, const char* file){

  double ** g_sequential = m_functionofinterest[ds]->get_g_sequential();
  ofstream outfile(file, ios::out);

  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }

  for(unsigned int i=0; i<m_num_of_intervals; i++){
    for(int j=0; j<m_length_grid; j++)
      outfile<<g_sequential[i][j]<<" ";
    outfile<<endl;
  }
  outfile.close();
}


void SMC_PP_MCMC::print_prob_sequential(int ds,const char* file){

  double ** prob_sequential = m_functionofinterest[ds]->get_prob_sequential();
  stringstream out;

  ofstream outfile(file, ios::out);
  outfile<<setiosflags(ios::fixed);
  outfile.precision(5);

  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }
  for(unsigned int i=0; i<m_num_of_intervals; i++){
    for(int j=0; j<m_length_grid; j++)
      outfile<<prob_sequential[i][j]<<" ";
    outfile<<endl;
  }
  outfile.close();
}


void SMC_PP_MCMC::print_prob(const char* file){
  ofstream outfile(file, ios::out);
  outfile<<setiosflags(ios::fixed);
  outfile.precision(5);
  
  if(!outfile){
    cerr<<outfile << " could not be opened"<<endl;
    return;
  }

  for(int ds=0; ds<m_num; ds++){
    long  double * prob = m_functionofinterest[ds]->get_prob();
    for(int i=0; i<m_length_grid; i++){
      outfile<<prob[i]<<" ";
    }
    outfile<<endl;
  }
  outfile.close();
}


void SMC_PP_MCMC::print_variable_sample_sizes(const char* file){
    ofstream outfile(file, ios::out);
    if(!outfile){
      cerr<<outfile << " could not be opened"<<endl;
      return;
    }

    for(int ds=0; ds<m_num; ds++){
      for(unsigned int i=0; i<m_num_of_intervals; i++)
	outfile<<m_sample_sizes[ds][i]<<" ";
      outfile<<endl;
    }

    outfile.close();
}

void SMC_PP_MCMC::set_look_back(bool lb){
  m_do_SMC_past=lb;

  if(!m_functionofinterest && m_do_SMC_past){
    for(int ds=0; ds<m_num; ds++){
      m_functionofinterest = new Function_of_Interest * [m_num];
      m_functionofinterest[ds] = new Function_of_Interest(m_length_grid,m_start,m_end,m_nu,1,0,m_calculate_intensity,1,0,m_num_of_intervals,0,0);
    }
  }
}
            
    
