#ifndef REJECTION_SAMPLING_HPP
#define RJECTION_SAMPLING_HPP

#include "particle.hpp"
#include "changepoint.hpp"
#include "probability_model.hpp"
#include "histogram_type.hpp"

/*allows for rejection sampling when there is only one changepoint*/

class rejection_sampling {
  
public:
  rejection_sampling(double m_start_time = 0, double m_cp_start = 0, double m_end_time = 1, long long int m_sample_size = 10, probability_model *m_pm = NULL,  double m_cp_prior = 1, int m_seed = 0, bool m_calculate_mean = 0, bool m_calculate_sample_histogram = false, unsigned int m_num_bins = 100);
  ~rejection_sampling();
  void write_1d_histogram_to_file(const string output_filename="1d_distribution.txt"){m_histogram->write_1d_histogram_to_file(output_filename);}
  void run_simulation();
  void sample_from_prior();
  Particle<changepoint>** get_sample() {return m_sample;}
  void write_frequency_counts_to_file(const char*);
  double m_acceptance_rate;


private:
  
  double m_start_time;
  double m_cp_start;
  double m_end_time;
  long long int m_sample_size;
  probability_model *m_pm;
  double m_cp_prior;
  int m_seed;
  const gsl_rng_type *m_r_type;
  gsl_rng *m_r;
  double m_zero_cp_likelihood;
  Particle<changepoint> **m_sample;
  double draw_from_prior();
  double calculate_mle(double, double);
  map<unsigned int,unsigned long long int> m_dimension_frequency_count;
  Histogram_Type<changepoint> *m_histogram;
  bool m_calculate_sample_histogram;
  int m_num_bins;
  bool m_calculate_mean;
  changepoint *m_end_of_int_changepoint;
  
};


#endif //REJECTION_SAMPLING_HPP
