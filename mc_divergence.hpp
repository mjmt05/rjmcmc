#ifndef MC_DIVERGENCE_HPP
#define MC_DIVERGENCE_HPP

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_cdf.h>

enum Divergence_Type { BIAS, CHI_SQUARE, FOI, EXTENT, ENTROPY, NONE };
enum Loss_Function { MINIMAX, AVERAGE };

class mc_divergence{

  
 public:
  mc_divergence(Divergence_Type=BIAS,Loss_Function=MINIMAX,unsigned int=0);
  virtual ~mc_divergence();
  void reset();
  void update_divergence(unsigned int);
  double get_divergence();
  const Divergence_Type get_divergence_type(){ return m_divergence_type; }
  const Loss_Function get_loss_function(){ return m_loss_fn; }
  void set_var(double var){ m_var = var; }
  void set_initial_number_of_bins();
  static void set_delta(double);
  static double calculate_chi_squared_criteria( unsigned int );
  static double phi(unsigned int);
  static double phi1(unsigned int);//first forward difference of phi()
  static double phi2(unsigned int);//second forward difference of phi()
  static double log1(unsigned int);//first forward difference of log()
  static void create_phi_lookup( unsigned int, Loss_Function=MINIMAX);
  static void delete_phi_lookup();
  static void create_chi_lookup( unsigned int );
  static void delete_chi_lookup();
  static void create_log_lookup( unsigned int );
  static void delete_log_lookup();
  static void create_lookup_arrays(unsigned int length, Divergence_Type=BIAS,Loss_Function=MINIMAX);
  static void delete_lookup_arrays(Divergence_Type=BIAS);

 protected:

  double m_divergence;
  double m_existing_bins_divergence;
  double m_new_bin_divergence_change;
  double m_new_bin_prob;
  unsigned int m_sample_size, m_initial_sample_size;
  Divergence_Type m_divergence_type;
  Loss_Function m_loss_fn;
  double m_var;
  unsigned long long int m_non_empty_bins,  m_initial_non_empty_bins;
  unsigned long long int m_waiting_time, m_last_waiting_time;//time since last new bin discovery
  static double* phi_lookup;
  static double* phi1_lookup;
  static double* phi2_lookup;
  static double* chi_lookup;
  static double* log_lookup;
  static double m_delta;
  static unsigned int length_phi_lookup;
  static unsigned int length_chi_lookup;
  static unsigned int length_log_lookup;
  static double m_phi_1;

};




#endif
