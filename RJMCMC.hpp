#ifndef RJMCMC_HPP
#define RJMCMC_HPP


#include "particle.hpp"
#include "histogram_type.hpp"
#include "mc_divergence.hpp"
using namespace std;
#define LOG_TWO log(2.0)
#define LOG_THREE log(3.0)
#define ONE_THIRD 1.0/3.0;
#include "Data.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <time.h>





template< class T>
  class rj{

 public:
  rj(double = 0, double =1, long long int = 10, int = 1000, long long int =1, long long int=100, bool=0,Particle<T> * = NULL, int seed = 0, bool = false);
  virtual ~rj();

  virtual void initiate_sample(Particle<T> *);
  virtual T* generate_new_parameter() const = 0;
  virtual T* move_parameter(unsigned int) const {return NULL;}  
  virtual T* copy_parameter(T *)const=0;
  virtual double log_likelihood_ratio_birth(T *, int) = 0;
  virtual double log_likelihood_ratio_death(unsigned int) = 0; 
  virtual double log_likelihood_ratio_move(T *, unsigned int) = 0;
  virtual double log_likelihood_ratio_move_parameter(int) = 0;
  virtual double log_prior_ratio_birth(T *) const = 0;
  virtual double log_prior_ratio_death(unsigned int) const = 0;
  virtual double log_prior_ratio_move(T *, unsigned int) const = 0;
  virtual double log_prior_ratio_move_parameter() const = 0;
  virtual double log_proposal_ratio_birth(T*) const = 0;
  virtual double log_proposal_ratio_death(int) const = 0;
  virtual double log_proposal_ratio_move() const = 0;
  virtual double log_proposal_ratio_move_parameter() const = 0;
  double birth_proposal(unsigned int) const;
  double death_proposal(unsigned int) const;
  double move_proposal(unsigned int) const;
  virtual void update_likelihood_birth(T*,int,double) = 0;
  virtual void update_likelihood_death(unsigned int, double) = 0;
  virtual void update_likelihood_move(T*, int, double) = 0;
  virtual void update_likelihood_move_parameter(int, double) = 0;
  virtual Particle<T> * copy_particle(Particle<T> *);
  virtual void destroy_sample();
  virtual void update_log_posterior();
  virtual Particle<T>** get_sample() const {return m_sample;}
  virtual Particle<T>* get_current_particle() const {return m_current_particle;}
  void delete_current_particle(){delete m_current_particle;}
  virtual double get_end_time() const {return m_end_time;}
  virtual void set_initial_sample(Particle<T> * initial) {m_initial_sample = initial;}
  virtual void calculate_function_of_interest() = 0;
  void open_sample_stream();
  void print_sample();
  void open_size_of_sample_stream();
  void print_size_of_sample();
  void open_logposterior_of_sample_stream();
  void print_logposterior_of_sample();
  void write_frequency_counts_to_file(const char* ="dimension_distribution.txt");
  void write_1d_histogram_to_file(const string output_filename="1d_distribution.txt"){m_histogram->write_1d_histogram_to_file(output_filename);}
  void write_marginal_histogram_to_file(const string output_filename="marginal_distribution.txt",unsigned int dimension=1,unsigned int margin=0){m_histogram->write_marginal_histogram_to_file(output_filename,dimension,margin);}
  void write_primary_function_of_interest_to_file(const string ="foi.txt");
  map<vector<unsigned int>,long long unsigned int>* get_bin_counts() {if (m_histogram) return m_histogram->get_bin_counts(); return NULL;}
  Histogram* get_histogram(){return m_histogram;}
  double get_autocorrelation_prob(){ return m_histogram->get_autocorrelation_prob();}
  unsigned int get_MAP_dimension() const {return m_MAP_dimension;};
  Particle<T> * get_MAP(unsigned int k) {return m_MAPs[k];};
  Particle<T> * get_MAP_dimension_MAP() {return m_MAPs[m_MAP_dimension];};
  virtual void locally_perfect_MAP(Particle<T>* = NULL, bool = false ) = 0;
  unsigned long long int get_size_sample()const{return m_size_of_sample;} 
  void set_continue_loop(bool loop){m_continue_loop=loop;}
  void runsimulation();
  void start_calculating_mc_divergence(Divergence_Type=BIAS,Loss_Function=MINIMAX,unsigned int=50,bool=false,unsigned int=0,bool=false,unsigned int=0);
  void update_dimension_distribution();
  virtual void update_histograms();
  virtual void update_primary_function_of_interest(){;}
  double get_divergence();
  void set_mc_divergence_chi_delta(double delta){ mc_divergence::set_delta(delta); }
  void end_divergence_burn_in(){ m_histogram->set_divergence_initial_number_of_bins();}
  void set_initial_iterations(long long int ii){m_initial_iterations=ii;}
  void set_function_criteria(unsigned int);
  void set_sample_filename(string s){m_sample_filename = s;}
  void set_sample_dimensions_filename(string s){m_sample_dimensions_filename = s;}
  void set_sample_logposterior_filename(string s){m_sample_logposterior_filename = s;}
  void start_printing_sample(bool v = 1, bool d = 1, bool p = 1);
  void stop_printing_sample();
  void start_storing_sample(){m_storing_sample=true; if(!m_sample) m_sample = new Particle<T>* [m_iterations];}
  void stop_storing_sample(){m_storing_sample=false;}
  void print_current_sample();
  virtual void print_acceptance_rates();
  void do_hill_climbing(){m_hill_climbing=true;}
  virtual void update_new_parameter_with_position( T*, unsigned int ) {}
  void calculate_sample_histogram( bool mv=false, unsigned int num_bins=50, bool bounded=true, bool weighted_histogram=false, bool calculate_divergence=false, Divergence_Type divergence_type=BIAS, Loss_Function loss_fn=MINIMAX, unsigned int look_up_length=0, bool estaimte_cor=false, unsigned int max_dim=0 );
  void calculate_current_bin();
  void no_1d_histogram(){ m_calculate_1d_sample_histogram = false; }
  virtual double get_component_value( T* ){ return 0;}
  virtual void set_histogram_component_value_function()=0;
  void do_importance_sampling(){ m_importance_sampling = true; }
  virtual void calculate_importance_weight(){}
  double get_sum_thinned_importance_weights(){ return m_sum_thinned_importance_weights; }
  void non_conjugate(){m_conjugate=0;}
  void constrained_model(){m_constraint=true;}
  virtual bool draw_means_from_posterior( bool test_constraint= true ){return true;}

 protected:
 
  double m_start_time;
  double m_end_time;
  long long int m_iterations;
  unsigned int m_max_theta;
  long long int m_thinning;
  long long int m_burnin;
  Particle<T> * m_initial_sample;
  unsigned int m_length_grid; 
  bool m_calculate_div;
  int m_seed;
  unsigned long long int m_size_of_sample;
  Particle<T> * m_current_particle;
  unsigned int m_k;//dimension of m_current_particle
  double m_log_likelihood_ratio;
  double m_log_proposal_ratio;
  double m_log_prior_ratio;
  double m_log_posterior_ratio;
  double m_b; //probability of doing a birth
  double m_d; //probability of doing a death
  double m_m; //probability of doing a move
  map<unsigned int,Particle<T> *> m_MAPs;//m_MAPs[k] gives the MAP model of dimension k.
  map<unsigned int,unsigned long long int> m_dimension_frequency_count;//stores frequency counts for the parameter dimensions visited during MCMC.
  map<unsigned int,long double> m_dimension_frequency_weight;//a weighted version, for importance sampling
  Histogram_Type<T> *m_histogram;
  unsigned int m_MAP_dimension;//the parameter dimension most frequently visited.
  bool m_calculate_sample_histogram;//determines whether or not to use dimension_frequency_count and map_dimension;
  bool m_calculate_mv_sample_histogram;//(only if m_calculate_sample_histogram) Determines whether or not to calcluate transdimensional histogram.
  bool m_calculate_1d_sample_histogram;//(only if m_calculate_sample_histogram) Determines whether or not to calcluate 1d histogram.
  bool m_calculate_function;
  double *m_mean_function_of_interest;
  double *m_mean_sq_function_of_interest;
  double *m_var_function_of_interest;
  double m_mean_var_function_of_interest;
  bool m_delete_data;
  unsigned long long int m_birth_attempt;
  unsigned long long int m_death_attempt;
  unsigned long long int m_move_attempt;
  unsigned long long int m_move_parameter_attempt;
  unsigned long long int m_birth_accept;
  unsigned long long int m_death_accept;
  unsigned long long int m_move_accept; 
  unsigned long long int m_move_parameter_accept;
  bool m_accept;
  bool m_accept_between_thinning;
  unsigned long long int m_acceptances_between_thinning;
  bool m_hill_climbing;//if true, accept iff logposterior increases
  bool m_calculating_divergence;
  bool m_importance_sampling;
  long double m_current_log_importance_weight, m_current_importance_weight;
  long double m_sum_importance_weights, m_sum_thinned_importance_weights;
  long long int m_initial_iterations;
  const gsl_rng_type * r_type;
  gsl_rng * r;
  bool m_continue_loop;
  long long int m_iters;
  Particle<T> ** m_sample;
  string m_sample_filename;
  string m_sample_dimensions_filename;
  string m_sample_logposterior_filename;
  ofstream m_sample_stream;
  ofstream m_sample_dimensions_stream;
  ofstream m_sample_logposterior_stream;
  bool m_storing_sample;
  double m_start_cps;
  bool m_printing_sample;
  bool m_printing_sample_value;//only valid if m_printing_sample == true
  bool m_printing_sample_dimension;//only valid if m_printing_sample == true
  bool m_printing_sample_posterior;//only valid if m_printing_sample == true
  bool m_conjugate;
  bool m_constraint;
  unsigned long long int m_num_constrained_particles;
  void rj_construct();
};



template<class T>
rj<T>::rj(double begin, double end, long long int its, int max, long long int thin, long long int burnin, bool calcKL,Particle<T> * initial,int se, bool store_sample)
:m_start_time(begin),m_end_time(end),m_iterations(its),m_max_theta(max),m_thinning(thin),m_burnin(burnin),m_initial_sample(initial),m_calculate_div(calcKL),m_seed(se),m_storing_sample(store_sample),m_start_cps(begin)
  {
    set_initial_sample(initial);
    rj_construct();
  }



template<class T>
  rj<T>::~rj()                                                                                
{
    
  if(m_sample!=NULL){
    destroy_sample();                   
    delete [] m_sample;
  }

  if (m_current_particle){
    delete m_current_particle;
  }

  for(typename std::map<unsigned int,Particle<T> *>::iterator iter = m_MAPs.begin(); iter != m_MAPs.end(); ++iter){
    delete iter->second;
  }

  m_MAPs.erase(m_MAPs.begin(),m_MAPs.end());

 
   gsl_rng_free(r);
  
  if(m_histogram)
    delete m_histogram;

  if(m_mean_function_of_interest)
    delete [] m_mean_function_of_interest;
  if(m_mean_sq_function_of_interest)
    delete [] m_mean_sq_function_of_interest;
  if(m_var_function_of_interest)
    delete [] m_var_function_of_interest;

}

template<class T>
  void rj<T>::rj_construct(){
  m_conjugate=true;
  m_constraint=false;
  m_num_constrained_particles=0;
  m_calculating_divergence = false;
  m_calculate_sample_histogram = false;
  m_calculate_mv_sample_histogram = false;
  m_calculate_1d_sample_histogram = true;

  m_continue_loop=0;
 
  m_sample = m_storing_sample ? new Particle<T>* [m_iterations] : NULL;
  m_current_particle = NULL;
  m_MAP_dimension = 0;
  m_histogram = NULL;

  if(m_calculate_div){
    m_initial_iterations=500;
  }

  m_mean_function_of_interest=NULL;
  m_mean_sq_function_of_interest=NULL;
  m_var_function_of_interest=NULL;
  m_length_grid=1;  

  m_size_of_sample=m_iterations;
   
  int seed=m_seed;
  gsl_rng_env_setup();
  r_type=gsl_rng_default;
  r = gsl_rng_alloc(r_type);
  gsl_rng_set(r,seed);

  m_sample_filename = "sampleRJ.txt";   
  m_sample_dimensions_filename = "sizesampleRJ.txt";
  m_sample_logposterior_filename = "logposterior_sampleRJ.txt";
  m_printing_sample = false;
  m_printing_sample_value = false;
  m_printing_sample_dimension = false;
  m_printing_sample_posterior = false;
  m_hill_climbing = false;
  m_importance_sampling = false;
  m_current_log_importance_weight = 0;
  m_current_importance_weight = 1;
  m_sum_importance_weights = m_iterations*m_thinning;
  m_sum_thinned_importance_weights = m_iterations;
}

template<class T>
void rj<T>::set_function_criteria(unsigned int grid){

  m_length_grid=grid;
  m_mean_function_of_interest =  new double [m_length_grid];
  m_mean_sq_function_of_interest = new double [m_length_grid];
  m_var_function_of_interest =  new double [m_length_grid];
  for(unsigned int i=0; i<m_length_grid; i++)
    m_mean_function_of_interest[i] = m_mean_sq_function_of_interest[i] = m_var_function_of_interest[i] = 0;

}



template<class T>
  void  rj<T>::runsimulation()
{
  
  bool first=0;

  //if this is the first time runsimulation is called initialize variables, allows you to go in an out of the loop used for variable number of particles
  if(!m_continue_loop){
    m_birth_attempt = 0;
    m_death_attempt = 0;
    m_move_attempt = 0;
    m_move_parameter_attempt=0;
    m_birth_accept = 0;
    m_death_accept = 0;
    m_move_accept = 0;
    m_move_parameter_accept=0;
    m_continue_loop=1;

    initiate_sample(m_initial_sample);
    m_iters = -m_burnin;
    first=1;
   
  }
  long long int n_iters;
  double u1,u2;
  bool move_type=0;
  T* new_value=NULL;
  T* new_value_alt=NULL;
  unsigned int new_position=0;
  int index_theta_delete=0;
  int index_theta_move=0;
  bool gsl_b=0;

  m_k = m_current_particle->get_dim_theta();
  m_accept = true;
  m_accept_between_thinning = false;
  m_acceptances_between_thinning = 0;
  //for (long long int i= -m_burnin; i<m_iterations*m_thinning; i++){
  while(m_continue_loop){

  
    n_iters=m_iters/m_thinning;
    
    if( m_accept ){
      m_b = birth_proposal(m_k);
      m_d = death_proposal(m_k);
      m_m = move_proposal(m_k);
    }

    m_log_proposal_ratio = 0;
   
    u1=gsl_ran_flat(r,0,1);
    /*birth*/
    if(u1<=m_b){
      new_value = generate_new_parameter();
      new_position = m_current_particle->find_position(new_value);
      if(new_position<=m_k){
	if(!m_conjugate){
	  new_value_alt=copy_parameter(m_current_particle->get_theta_component(new_position-1));
	}
	update_new_parameter_with_position( new_value, new_position );
	m_log_likelihood_ratio = log_likelihood_ratio_birth(new_value, new_position);
	m_log_prior_ratio = log_prior_ratio_birth(new_value);
	if(!m_k){
	  if(m_conjugate)
	    m_log_proposal_ratio = (m_max_theta == 1) ? -LOG_TWO:-LOG_THREE;
	  else
	    m_log_proposal_ratio = (m_max_theta == 1) ? 0: LOG_TWO-LOG_THREE;
	}
	else if(m_k == m_max_theta-1)
	  m_log_proposal_ratio = LOG_THREE-LOG_TWO;
	m_log_proposal_ratio += log_proposal_ratio_birth(new_value);
      }else{
	/* was having problems where gsl was generating the same double this ignores that iteration*/
	m_log_likelihood_ratio=-1;
	m_log_prior_ratio=-1;
	m_log_proposal_ratio=-1;
	gsl_b=1;
      }
      if(m_iters>0)
	m_birth_attempt++;

    

    }
    /*death*/
    else if (u1<=m_b+m_d && u1>m_b){
      index_theta_delete = static_cast<int>( gsl_ran_flat( r, 0, m_k ) );
      if(!m_conjugate){
	
	new_value=copy_parameter(m_current_particle->get_theta_component(index_theta_delete-1));
      }
      m_log_likelihood_ratio = log_likelihood_ratio_death(index_theta_delete);
      m_log_prior_ratio = log_prior_ratio_death(index_theta_delete);
      if(m_k == 1){
	if(m_conjugate)
	  m_log_proposal_ratio = (m_max_theta == 1) ? LOG_TWO:LOG_THREE;
	else
	  m_log_proposal_ratio = (m_max_theta == 1) ? 0:LOG_THREE-LOG_TWO;
      }
      else if(m_k == m_max_theta)
	m_log_proposal_ratio = LOG_TWO-LOG_THREE;
      m_log_proposal_ratio += log_proposal_ratio_death(index_theta_delete);
      if(m_iters>0)
	m_death_attempt++;
    }
    /*move*/
    else {
      if(!m_conjugate && (!m_k || gsl_ran_flat(r,0,1)<0.5)){
	index_theta_move = gsl_rng_uniform_int (r,m_k+1)-1;
	if(m_iters>0)
	  m_move_parameter_attempt++;


	new_value=copy_parameter(m_current_particle->get_theta_component(index_theta_move));

	m_log_likelihood_ratio = log_likelihood_ratio_move_parameter(index_theta_move);

	m_log_prior_ratio=log_prior_ratio_move_parameter();
	m_log_proposal_ratio=log_proposal_ratio_move_parameter();
	move_type=1;
      }else{

	index_theta_move = static_cast<unsigned int>( gsl_ran_flat( r, 0, m_k ) );
	new_value =  move_parameter(index_theta_move);
      	if(!m_conjugate)
	  new_value_alt=copy_parameter(m_current_particle->get_theta_component(index_theta_move-1));
	m_log_likelihood_ratio = log_likelihood_ratio_move(new_value, index_theta_move);
	m_log_prior_ratio = log_prior_ratio_move(new_value, index_theta_move);
	m_log_proposal_ratio += log_proposal_ratio_move();
	if(m_iters>0)
	  m_move_attempt++;

      }
    }

    m_log_posterior_ratio = m_log_likelihood_ratio + m_log_prior_ratio;

    if(m_iters>0 && m_hill_climbing){
      m_accept = m_log_posterior_ratio > 0;
    }
    else{
      if (m_log_posterior_ratio + m_log_proposal_ratio >= 0 ){
	m_accept = 1;
      }else{        
	if(gsl_b){
	  u2=0;
	  gsl_b=0;
	}else{
	  u2=log( gsl_ran_flat( r, 0, 1 ) );}
	  
	if (u2<m_log_posterior_ratio + m_log_proposal_ratio){
	  m_accept = 1;
	} else{
	  m_accept = 0;
	}
      }
    }
    if (m_accept) {
      m_accept_between_thinning = true;
      /*birth*/
      if (u1<=m_b){
	update_likelihood_birth(new_value,new_position,m_log_likelihood_ratio);
	m_current_particle -> add_component(new_value,new_position);
	delete new_value_alt;
	new_value_alt=NULL;
	if(m_iters>0)
	  m_birth_accept++;
      }
      /*death*/
      else if (u1<=m_b+m_d && u1>m_b){
	update_likelihood_death(index_theta_delete,m_log_likelihood_ratio);                
	m_current_particle -> delete_component(index_theta_delete);

	if(!m_conjugate){
	  delete new_value;
	  new_value=NULL;
	}

	if(m_iters>0)
	  m_death_accept++;
      }
      /*move*/
      else{

	if(move_type){
	  update_likelihood_move_parameter(index_theta_move,m_log_likelihood_ratio);
	  delete new_value;
	  new_value=NULL;
	  move_type=0;
	  if(m_iters>0)
	    m_move_parameter_accept++;
        
	}else{
	  update_likelihood_move(new_value,index_theta_move,m_log_likelihood_ratio);
	  m_current_particle -> change_component(new_value, index_theta_move);
	  delete new_value_alt;
	  new_value_alt=NULL;
	  if(m_iters>0)
	    m_move_accept++;
	}
      }

      m_k = m_current_particle->get_dim_theta();

      update_log_posterior();
      if(m_iters>0 && ((m_iters % m_thinning) == 0)){
	if(m_accept_between_thinning){
	  m_acceptances_between_thinning++;
	}
	m_accept_between_thinning = false;
	if(m_storing_sample)
	  m_sample[n_iters-1]=copy_particle(m_current_particle);
	else if(m_sample)
	  m_sample[n_iters-1]=m_current_particle;
	if(m_printing_sample)
	  print_current_sample();
      }
    }else{
      
      if (u1<m_b|| (u1>m_b+m_d && !move_type)){
	delete new_value;
      }

      if(!m_conjugate){
	if(u1<m_b){
	  m_current_particle->change_component(new_value_alt,new_position-1);
	  new_value_alt=NULL;
	}else if(u1<=m_b+m_d && u1>m_b){
	  m_current_particle->change_component(new_value,index_theta_delete-1);
	  new_value=NULL;
	}else{
	  if(move_type){
	    m_current_particle->change_component(new_value,index_theta_move);
	    new_value=NULL;
	    move_type=0;
	  }else{
	    m_current_particle->change_component(new_value_alt,index_theta_move-1);
	    new_value_alt=NULL;
	  }
	}
      }

              
      if(m_iters==m_thinning){
	if(m_storing_sample){
	  m_sample[0]= copy_particle(m_current_particle);
	}
	else if(m_sample){
	  m_sample[0]=m_current_particle;
	}
	if(m_printing_sample)
	  print_current_sample();
      }

      if(m_iters>m_thinning  && ((m_iters % m_thinning) == 0)){
	
	if(m_sample){
	  if (m_sample[n_iters -2] == m_current_particle){
	    m_sample[n_iters-1] = m_sample[n_iters -2];
	  }
	  else{
	    m_sample[n_iters-1]=copy_particle(m_current_particle);
	  }
	}
	if(m_printing_sample)
	  print_current_sample();
      }
    }

    bool particle_okay=true;
    if(m_constraint && m_iters>0 && m_iters%m_thinning==0){
      particle_okay=draw_means_from_posterior();
      if(particle_okay)
	m_num_constrained_particles++;
    }
    
    if(m_importance_sampling && m_iters>0){
      if(m_iters == 1 || m_accept ){
	if(particle_okay)
	  calculate_importance_weight();
	else
	  m_current_log_importance_weight=-DBL_MAX;
	cout << m_current_log_importance_weight << endl;
        m_current_particle->set_weight(m_current_log_importance_weight);
        m_current_importance_weight = exp(m_current_log_importance_weight);
      }
      m_sum_importance_weights += m_current_importance_weight;
      if(!(m_iters%m_thinning))
        m_sum_thinned_importance_weights += m_current_importance_weight;
    }
        
    if(m_iters>0 && particle_okay)
      update_dimension_distribution();

    if(m_iters>0 && m_iters%m_thinning==0 && particle_okay)
    {
      if(m_mean_function_of_interest){
	update_primary_function_of_interest();
	if(m_calculating_divergence){
	  m_mean_var_function_of_interest = 0;
	  for(unsigned int i = 0; i < m_length_grid; i++ ){
	    double var_mean_i =  (m_mean_sq_function_of_interest[i] - m_mean_function_of_interest[i]*(m_mean_function_of_interest[i]/n_iters))/(n_iters-1);
	    double mean_var_i = m_var_function_of_interest[i]/n_iters;
	    m_mean_var_function_of_interest += var_mean_i + mean_var_i;
	  }
	  m_mean_var_function_of_interest /= m_length_grid;
	  m_histogram->set_divergence_var(m_mean_var_function_of_interest);
	}
      }
      if(m_calculate_function && !m_storing_sample){
	calculate_function_of_interest();
      }
    }

    if(m_calculate_sample_histogram && m_iters>0 && particle_okay)
      update_histograms();


    if(m_calculate_div && m_iters>0 && m_iters%m_thinning==0){

      if(n_iters%(m_initial_iterations)==0){
	m_size_of_sample=n_iters;
	m_continue_loop=0;	
	if(first){
	  m_initial_iterations=1;
	  first=0;
	}
      }
    }

   
    m_iters++;

    if(m_iters>m_iterations*m_thinning){
      m_continue_loop=0;
    }
  }
}


template<class T>
void rj<T>::initiate_sample(Particle<T> * ptr2particle){
  if (ptr2particle==NULL){
    m_current_particle = new Particle<T>(0,NULL);
  }
  else{
    m_current_particle = ptr2particle;
  }
}


template<class T>
Particle<T> * rj<T>::copy_particle(Particle<T> * ptr2particle){
  Particle<T> * new_particle;
  new_particle = new Particle<T>(ptr2particle,NULL);
  return(new_particle);
}

template<class T>
double rj<T>::birth_proposal(unsigned int k) const{
  double b;
    
  if (k==0){
    b=m_conjugate?1:.5;
  }
  else if(k==m_max_theta){
    b=0;
  }
  else{
    b=ONE_THIRD;
  }

  return (b);
}

template<class T>
  double rj<T>::death_proposal(unsigned int k) const {
  double d;

  if (k==0){
    d=0;
  }
  else if(k==m_max_theta){
    d=0.5;
  }
  else{
    d=ONE_THIRD;
  }

  return (d);
}

template<class T>
double rj<T>::move_proposal(unsigned int k) const{ 

  double m;

  if(k==0){
    m=m_conjugate?0:.5;
  }
  else if(k==m_max_theta){
    m=0.5;
  }
  else{
    m=ONE_THIRD;
  }

  return (m);
}

    
template<class T>    
  void rj<T>::update_log_posterior(){
  m_current_particle->set_log_posterior(m_current_particle->get_log_posterior()+m_log_posterior_ratio);
}


template<class T>
void rj<T>::destroy_sample(){
 
  for (long long int i=m_size_of_sample-1; i>0; i--){
 
    if (m_sample[i]!=m_sample[i-1]){    
      delete m_sample[i];
    }
  }


  if( m_sample[0] != m_current_particle )
    delete m_sample[0];

}



template<class T>
void rj<T>::open_sample_stream(){
  m_sample_stream.open(m_sample_filename.c_str(),ios::out);
  if(!m_sample_stream){
    cerr<<"Sample file "<<m_sample_filename.c_str()<<" could not be opened"<<endl;
    return;
    //exit(1);
  }
  m_sample_stream<<setiosflags(ios::fixed);
  m_sample_stream.precision(5);
}

template<class T>
void rj<T>::print_sample(){

  open_sample_stream();
  for(long long int i=0; i<m_size_of_sample; i++){
      m_sample_stream<<*m_sample[i]<<endl;
  }
  m_sample_stream.close();
}

template<class T>
void rj<T>::open_size_of_sample_stream(){
  
  m_sample_dimensions_stream.open(m_sample_dimensions_filename.c_str(), ios::out);
  
  if(!m_sample_dimensions_stream){
    cerr << "Dimensions file " << m_sample_dimensions_filename.c_str() << "could not be opened" << endl;
    return;
    //exit(1);
  }
}

template<class T>
  void rj<T>::print_size_of_sample(){

  open_size_of_sample_stream();
  for(long long int i=0; i<m_size_of_sample; i++){
    m_sample_dimensions_stream<<m_sample[i]->get_dim_theta()<<' ';
  }

  m_sample_dimensions_stream.close();
}

template<class T>
  void rj<T>::open_logposterior_of_sample_stream(){
  
  m_sample_logposterior_stream.open(m_sample_logposterior_filename.c_str(), ios::out);
  
  if(!m_sample_logposterior_stream){
    cerr<<"Log posterior file " << m_sample_logposterior_filename.c_str() << " could not be opened" << endl;
    return;
    //exit(1);
  }
}

template<class T>
  void rj<T>::print_logposterior_of_sample(){

  open_logposterior_of_sample_stream();
  for(long long int i=0; i<m_size_of_sample; i++){
    m_sample_logposterior_stream<<m_sample[i]->get_log_posterior()<<' ';
  }

  m_sample_logposterior_stream.close();
}

template<class T>
void rj<T>::start_printing_sample(bool value, bool dimension, bool posterior){
  m_printing_sample = true;
  m_printing_sample_value = value;
  m_printing_sample_dimension = dimension;
  m_printing_sample_posterior = posterior;
  
  if(m_printing_sample_value)
    open_sample_stream();
  if(m_printing_sample_dimension)
    open_size_of_sample_stream();
  if(m_printing_sample_posterior)
    open_logposterior_of_sample_stream();
}

template<class T>
  void rj<T>::stop_printing_sample(){
  m_printing_sample = false;
  if(m_sample_stream.is_open())
    m_sample_stream.close();
  if(m_sample_dimensions_stream.is_open())
    m_sample_dimensions_stream.close();
  if(m_sample_logposterior_stream.is_open())
    m_sample_logposterior_stream.close();
}

template<class T>
  void rj<T>::print_current_sample(){
  if( m_sample_stream.is_open() )
    m_sample_stream<<*m_current_particle;
  if( m_sample_dimensions_stream.is_open() )
    m_sample_dimensions_stream<<m_k<<" ";
  if( m_sample_logposterior_stream.is_open() )
    m_sample_logposterior_stream<<m_current_particle->get_log_posterior()<<" ";
}

template<class T>
  void rj<T>::print_acceptance_rates(){
  cout << "\nAcceptance Rates:" << endl;
  if( m_birth_attempt )
    cout << "\tBirth: " << m_birth_accept/(double)m_birth_attempt<<endl;
  if( m_death_attempt )
    cout << "\tDeath: " << m_death_accept/(double)m_death_attempt<< endl;
  if( m_move_attempt )
    cout << "\tMove:  " << m_move_accept/(double)m_move_attempt <<endl;
  cout << "\tOverall acceptance for thinned samples: " << m_acceptances_between_thinning/(double)m_size_of_sample << endl;
  if(!m_conjugate){
    if(m_move_parameter_attempt)
      cout << "\tMove Parameter:  " << m_move_parameter_accept/(double)m_move_parameter_attempt <<endl;
  }
  if(m_constraint)
    cout << "\tP(Constraint): " <<  m_num_constrained_particles/(double)m_size_of_sample << endl;
}

template<class T>
  void rj<T>::write_frequency_counts_to_file(const char* output_filename){
  ofstream OutputStream(output_filename, ios::out);
  if(!m_importance_sampling){
      OutputStream<<setiosflags(ios::fixed);
      OutputStream.precision((int)log10(m_iterations*m_thinning));
      map<unsigned int,unsigned long long int>::iterator iter = m_dimension_frequency_count.begin();
      while(iter != m_dimension_frequency_count.end()){
          if(iter->second)
              OutputStream << iter->first << "\t" << (double)iter->second/m_iters << endl;
          ++iter;
      }
  }else{
      map<unsigned int,long double>::iterator iter = m_dimension_frequency_weight.begin();
      while(iter != m_dimension_frequency_weight.end()){
          if(iter->second>0)
              OutputStream << iter->first << "\t" << (double)iter->second/m_sum_importance_weights << endl;
          ++iter;
      }
  }
  OutputStream.close();
}

template<class T>
void rj<T>::calculate_sample_histogram( bool mv, unsigned int num_bins, bool bounded, bool weighted_histogram, bool calculate_divergence, Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int look_up_length, bool estimate_cor, unsigned int max_dim ){ 
  m_calculate_sample_histogram = true;
  if(mv)
    m_calculate_mv_sample_histogram = true;

  double bin_width = (m_end_time-m_start_cps)/num_bins;
  
  if( !m_histogram ){
    if(!calculate_divergence){
      m_histogram = new Histogram_Type<T>( m_start_cps, m_end_time, num_bins, bin_width, bounded, weighted_histogram, m_calculate_1d_sample_histogram, max_dim );
    }
    else{
      m_histogram = new Histogram_Type<T>( m_start_cps, m_end_time, num_bins, bin_width, bounded, weighted_histogram, m_calculate_1d_sample_histogram, max_dim, NULL, true, divergence_type, loss_fn, look_up_length, estimate_cor );
      m_histogram->track_entropy();
//      m_histogram->use_single_component();
    }
    set_histogram_component_value_function();
  }
}

template<class T>
void rj<T>::start_calculating_mc_divergence(Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int num_bins, bool bounded, unsigned int look_up_length, bool estimate_cor, unsigned int max_dim){
  m_calculating_divergence = true;
  if(!m_calculate_sample_histogram)
    calculate_sample_histogram(true,num_bins,bounded,false,true,divergence_type,loss_fn,look_up_length,estimate_cor,max_dim);
  else
    m_calculate_mv_sample_histogram = true;
}

template<class T>
  void rj<T>::update_histograms(){
  bool end_of_batch_should_be_new = (m_thinning>1 && m_iters%m_thinning==0) || (m_thinning==1 && m_accept );
  bool end_of_batch = m_iters%m_thinning==0;
  if((m_iters==1)||(m_calculate_1d_sample_histogram && m_accept)||end_of_batch_should_be_new)
    m_histogram->calculate_bin( m_current_particle, m_k );
  else if(m_iters>1&&(m_thinning==1 && !m_accept))
    m_histogram->increment_num_bin_repeats();
  if(m_calculate_mv_sample_histogram && end_of_batch)
    m_histogram->increment_bin_counts(NULL,m_current_importance_weight);
  if(m_calculate_1d_sample_histogram)
    m_histogram->increment_1d_bin_counts(m_current_importance_weight);
}

template<class T>
  double rj<T>::get_divergence(){
  double divergence = m_histogram->get_divergence();
  return divergence;
}

template<class T>
  void rj<T>::write_primary_function_of_interest_to_file(const string output_filename){
  if(!m_mean_function_of_interest)
    return;
  unsigned int iterations = m_constraint? m_num_constrained_particles:m_iters/m_thinning;
  ofstream OutputStream(output_filename.c_str(), ios::out);
  for(unsigned int i = 0; i < m_length_grid; i++ ){
    double var_mean_i = (m_mean_sq_function_of_interest[i] - m_mean_function_of_interest[i]*(m_mean_function_of_interest[i]/iterations))/(m_iters-1);
    double mean_var_i = m_var_function_of_interest[i]/iterations;
    OutputStream << m_mean_function_of_interest[i]/iterations << "\t" << var_mean_i + mean_var_i << endl;
  }
  OutputStream.close();
}

template<class T>
  void rj<T>::update_dimension_distribution(){
    unsigned int current_dimension = m_k;
    if(!m_importance_sampling){
      bool entered_new_dimension = m_dimension_frequency_count.find(current_dimension) == m_dimension_frequency_count.end();
      if(entered_new_dimension || (m_MAPs[current_dimension]->get_log_posterior()<m_current_particle->get_log_posterior())){
	if(!entered_new_dimension && m_MAPs[current_dimension])
	  delete m_MAPs[current_dimension];
	m_MAPs[current_dimension] = copy_particle(m_current_particle);
      }
      if(entered_new_dimension)
        m_dimension_frequency_count[current_dimension] = 1;
      else
        m_dimension_frequency_count[current_dimension]++;
      if(m_MAP_dimension != current_dimension && ((m_dimension_frequency_count.find(m_MAP_dimension)==m_dimension_frequency_count.end())|| (m_dimension_frequency_count[current_dimension]>m_dimension_frequency_count[m_MAP_dimension])))
	m_MAP_dimension = current_dimension;
    }
    else{
      bool entered_new_dimension = m_dimension_frequency_weight.find(current_dimension) == m_dimension_frequency_weight.end();
      if(entered_new_dimension || (m_MAPs[current_dimension]->get_log_posterior()+m_MAPs[current_dimension]->get_weight()<m_current_particle->get_log_posterior()+m_current_particle->get_weight())){
        if(!entered_new_dimension && m_MAPs[current_dimension])
	  delete m_MAPs[current_dimension];
	m_MAPs[current_dimension] = copy_particle(m_current_particle);
      }
      if(entered_new_dimension)
	m_dimension_frequency_weight[current_dimension] = m_current_importance_weight;
      else
	m_dimension_frequency_weight[current_dimension] += m_current_importance_weight;
      if(m_MAP_dimension != current_dimension && ((m_dimension_frequency_weight.find(m_MAP_dimension)==m_dimension_frequency_weight.end())|| (m_dimension_frequency_weight[current_dimension]>m_dimension_frequency_weight[m_MAP_dimension])))
	m_MAP_dimension = current_dimension;
    }
}


#endif
