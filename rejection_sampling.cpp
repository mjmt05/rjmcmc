#include "rejection_sampling.hpp"
#include "Data.hpp"
#include "changepoint.hpp"
#include "histogram_type.hpp"
#include <vector>
rejection_sampling::rejection_sampling(double start, double cp_start, double end, long long int sample_size, probability_model *pm, double cp_prior, int seed, bool calculate_mean, bool calculate_sample_histogram, unsigned int grid) {
  m_start_time = start;
  m_cp_start = cp_start;
  m_end_time = end;
  m_sample_size = sample_size;
  m_pm = pm;
  m_cp_prior = cp_prior;
  m_seed = seed;
  gsl_rng_env_setup();
  m_r_type = gsl_rng_default;
  m_r = gsl_rng_alloc(m_r_type);
  gsl_rng_set(m_r, seed);
  m_calculate_mean = calculate_mean;
  m_sample = NULL;
  m_dimension_frequency_count[0] = 0;
  m_dimension_frequency_count[1] = 0;
  m_histogram = NULL;
  m_calculate_sample_histogram = calculate_sample_histogram;
  m_num_bins = grid;
    
  if (m_calculate_sample_histogram == 1) {
    double bin_width = (m_end_time - m_cp_start) / (double) m_num_bins;
    bool calculate_1d_sample_histogram = true;
    bool weighted_histogram = false;
    bool bounded = true;
    int max_dim = 0;
    m_histogram = new Histogram_Type<changepoint>( m_cp_start, m_end_time, m_num_bins, bin_width, bounded, weighted_histogram, calculate_1d_sample_histogram, max_dim );
    m_histogram->set_value_function( &changepoint::getchangepoint );
  }

  m_end_of_int_changepoint = new changepoint(m_end_time, 0, 0, 0);
  m_pm->set_data_index(m_end_of_int_changepoint);
}

rejection_sampling::~rejection_sampling() {
  gsl_rng_free(m_r);
  if (m_sample) {
    for (long long int i = 0; i < m_sample_size; i++) {
      delete m_sample[i];
    }
    delete [] m_sample;
  }

  if (m_histogram) {
    delete m_histogram;
  }
  delete m_end_of_int_changepoint;
}


void rejection_sampling::sample_from_prior() {
  int ncps;
  double length_of_interval = m_end_time - m_cp_start;
  changepoint *cpintercept = NULL;
  changepoint *cp = NULL, *cp1 = NULL; 
  unsigned long long int data_index_start = m_pm->get_data()->find_data_index(m_cp_start);
//cout << current_likelihood << " " << data_index_end << " " << data_index_start << " " << m_start_time << " " << m_end_time;
  m_sample = new Particle<changepoint> *[m_sample_size];

  for (int i = 0; i < m_sample_size; i++) {
    ncps = gsl_ran_poisson(m_r, m_cp_prior * length_of_interval);
    cpintercept = new changepoint(m_cp_start, data_index_start, 0, 0);
    m_sample[i] = new Particle<changepoint>(0, NULL, cpintercept);
    if (ncps > 0) {
      vector<double> changepoints (ncps);
      for (int j = 0; j < ncps; j++) {
	changepoints[j] = gsl_ran_flat(m_r, m_cp_start, m_end_time);
      }
      sort(changepoints.begin(), changepoints.end());
      cp = new changepoint(changepoints[0], 0, 0, 0);
      m_pm->set_data_index(cp, 0, cpintercept);
      cpintercept->setlikelihood(m_pm->log_likelihood_interval(cpintercept, cp));
      cpintercept->setmeanvalue(m_pm->calculate_mean(cpintercept, cp));			 
      for (int j = 0; j < ncps - 1; j++) {
	cp1 = new changepoint(changepoints[j + 1], 0, 0, 0);
	m_pm->set_data_index(cp1, 0, cp);
	cp->setlikelihood(m_pm->log_likelihood_interval(cp, cp1));
	if (m_calculate_mean) {
	  cp->setmeanvalue(m_pm->calculate_mean(cp, cp1));
	}
	m_sample[i]->add_component(cp, j);
	cp = cp1;
      }
      cp->setlikelihood(m_pm->log_likelihood_interval(cp, m_end_of_int_changepoint));
      if (m_calculate_mean) {
	cp->setmeanvalue(m_pm->calculate_mean(cp, m_end_of_int_changepoint));
      }
      m_sample[i]->add_component(cp, ncps-1);
      changepoints.clear();
    } else {
      cpintercept->setlikelihood(m_pm->log_likelihood_interval(cpintercept, m_end_of_int_changepoint));
      cpintercept->setmeanvalue(m_pm->calculate_mean(cpintercept, m_end_of_int_changepoint));
    }	
  }    


}

void
rejection_sampling::run_simulation() {
  if (m_cp_start >= m_end_time) {
    cerr << "rejection_sampling: no valid changepoints" << endl;
    exit(1);
  }
  Data<double> *data = m_pm->get_data();
  unsigned long long int data_index_start = data->find_data_index(m_start_time);
  unsigned long long int data_index_end = m_end_of_int_changepoint->getdataindex();
  unsigned long long int data_index_temp = 0;
  double mle = calculate_mle(data_index_start, data_index_end);
  
  int accepted_samples = 0;
  double current_cp;
  double likelihood_left = 0, likelihood_right = 0;
  double log_posterior = 0;
  double u_rv;
  changepoint **cpvector = NULL; 
  changepoint *cpintercept = NULL;
  int n_cps = 0;
  m_sample =  new Particle<changepoint> *[m_sample_size];
  long long int attempts = 0;
  
  while (accepted_samples < m_sample_size) {
    attempts++;
    current_cp = draw_from_prior();
    if (current_cp > m_end_time) {
      likelihood_left = m_zero_cp_likelihood;
    } else {
      data_index_temp = data->find_data_index(current_cp, 0, data_index_start);
      likelihood_left = m_pm->log_likelihood_interval_with_count(m_start_time, current_cp, data_index_temp - data_index_start);
      likelihood_right = m_pm->log_likelihood_interval_with_count(current_cp, m_end_time, data_index_end - data_index_temp);
    }
    u_rv = log( gsl_ran_flat(m_r, 0, 1));

    if (likelihood_left + likelihood_right > mle) {
      cerr << "rejection_sampling: likelihood is greater than mle" << endl;
      exit(1);
    }

    if (u_rv < likelihood_left + likelihood_right - mle) {
      cpintercept = new changepoint(m_start_time, data_index_start, likelihood_left, 0);
     
      if (current_cp <= m_end_time) {
	n_cps = 1;
	cpvector = new changepoint *[1];
	cpvector[0] = new changepoint(current_cp, data_index_temp, likelihood_right, 0);
	m_dimension_frequency_count[1]++;
	if (m_calculate_mean) {
	  cpintercept->setmeanvalue(m_pm->calculate_mean(cpintercept, cpvector[0]));
	  cpintercept->setvarvalue(m_pm->get_var());
	  cpvector[0]->setmeanvalue(m_pm->calculate_mean(cpvector[0], m_end_of_int_changepoint));
	  cpvector[0]->setvarvalue(m_pm->get_var());
	}
      } else {
	m_dimension_frequency_count[0]++;
	if (m_calculate_mean) {
	  cpintercept->setmeanvalue(m_pm->calculate_mean(cpintercept, m_end_of_int_changepoint));
	  cpintercept->setvarvalue(m_pm->get_var());
	}
      }
      m_sample[accepted_samples] = new Particle<changepoint>(n_cps, cpvector, cpintercept);
      log_posterior = likelihood_left + likelihood_right;
      log_posterior += n_cps * log(m_cp_prior) - m_cp_prior * (m_end_time - m_cp_start);
      if (n_cps == 1) {
	log_posterior += m_cp_prior * (m_end_time - current_cp);
      }
      m_sample[accepted_samples]->set_log_posterior(log_posterior);
      if (m_calculate_sample_histogram) {
	m_histogram->calculate_bin(m_sample[accepted_samples], n_cps);
	m_histogram->increment_1d_bin_counts(1);
      }
      
      accepted_samples++;
      n_cps = 0;
      cpvector = NULL;
      cpintercept = NULL;

      
     //cout << "accept" <<endl;
    }
    likelihood_left = likelihood_right = 0;
  }
  m_acceptance_rate = (double) m_sample_size / attempts;
}

double 
rejection_sampling::draw_from_prior() {
  return gsl_ran_exponential (m_r, 1.0/m_cp_prior) + m_cp_start;
}

double 
rejection_sampling::calculate_mle(double data_index_start, double data_index_end) {
  Data<double> *data = m_pm->get_data();
  unsigned long long int data_index_1 = data->find_data_index(m_cp_start, 0, data_index_start);
  double current_likelihood;
  double current_data_point;
  double mle = 0;
  if (m_start_time != m_cp_start) {
    current_likelihood =  m_pm->log_likelihood_interval_with_count(m_start_time, m_cp_start, data_index_1 - data_index_start);
    current_likelihood += m_pm->log_likelihood_interval_with_count(m_cp_start, m_end_time, data_index_end - data_index_1);
    mle = current_likelihood;
  }

  for (int i = data_index_1; i < data_index_end; i++) {
    current_data_point = data->m_X[0][i];
    current_likelihood = m_pm->log_likelihood_interval_with_count(m_start_time, current_data_point, i - data_index_start);
    current_likelihood += m_pm->log_likelihood_interval_with_count(current_data_point, m_end_time, data_index_end - i);
    if (current_likelihood > mle || mle == 0) {
      mle = current_likelihood;
    }
    current_likelihood = m_pm->log_likelihood_interval_with_count(m_start_time, current_data_point + DBL_MIN, 
								  i + 1 - data_index_start);
    current_likelihood += m_pm->log_likelihood_interval_with_count(current_data_point + DBL_MIN, m_end_time, 
								   data_index_end - i - 1);
    if (current_likelihood > mle) {
      mle = current_likelihood;
    }
  }

  /* corresponds to no changepoint */
  current_likelihood = m_pm->log_likelihood_interval_with_count(m_start_time, m_end_time, data_index_end - data_index_start);
  //cout << current_likelihood << " " << data_index_end << " " << data_index_start << " " << m_start_time << " " << m_end_time;
  m_zero_cp_likelihood = current_likelihood;
  if (current_likelihood > mle) {
    mle = current_likelihood;
  }

  return mle;
  
}

void
rejection_sampling::write_frequency_counts_to_file(const char* output_filename) {
   ofstream OutputStream(output_filename, ios::out);
   OutputStream<<setiosflags(ios::fixed);
   //OutputStream.precision((int)log10());
   map<unsigned int,unsigned long long int>::iterator iter = m_dimension_frequency_count.begin();
   while(iter != m_dimension_frequency_count.end()){
     if(iter->second)
       OutputStream << iter->first << "\t" << (double)iter->second/m_sample_size << endl;
     ++iter;
   }
   OutputStream.close();

}
