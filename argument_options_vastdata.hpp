#ifndef ARGUMENT_OPTIONS_VASTDATA_HPP
#define ARGUMENT_OPTIONS_VASTDATA_HPP
#include <getopt.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include "mc_divergence.hpp"
using namespace std;

class ArgumentOptionsVast{
 public:
  ArgumentOptionsVast();
  ~ArgumentOptionsVast(){;};
  void parse(int argc, char *argv[]);
  double stringtodouble(char * parameter, char opt);
  long int stringtolong(char * parameter, char opt);
  void errorchecking(char * pend, char opt);
    
 
  /*mandatory arguments*/
  double m_start;
  double m_end; 
  string m_datafile; 
 

  /*optional arguments*/
  int m_num_intervals;
  unsigned long int m_particles;
  double m_cp_prior;
  double m_gamma_prior_1;
  double m_gamma_prior_2;
  int m_seed;
  bool m_calculate_filtering_mean;
  int m_grid;
  double m_ESS_threshold;
  bool m_print_ESS;
  unsigned int m_num_bins;
  Loss_Function m_loss_type;
  unsigned int m_min_iterations;
  bool m_fixed_sample_size;
  string m_sample_sizes;

  /*RJ paramters when sampling on the intervals over time*/
  int m_burnin;
  int m_thinning;
  double m_move_width;
  bool m_disallow_empty_intervals_between_cps;

 private:
  void usage(int status, char *);
  static struct option lopts[];

};

#endif
