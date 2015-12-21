#ifndef ARGUMENT_OPTIONS_HPP
#define ARGUMENT_OPTIONS_HPP
#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <string>
using namespace std;

class ArgumentOptions{
 public:
  ArgumentOptions();
  ~ArgumentOptions(){;};
  void parse(int argc, char *argv[]);
  double stringtodouble(char * parameter, char opt);
  long int stringtolong(char * parameter, char opt);
  void errorchecking(char * pend, char opt);
    
 
  /*mandatory arguments*/
  double m_start;
  double m_end; 
  string m_datafile; 
 

  /*optional arguments*/
  long int m_iterations;
  int m_burnin;
  int m_thinning;
  double m_move_width;
  double m_cp_prior;
  double m_gamma_prior_1;
  double m_gamma_prior_2;
  int m_seed;
  bool m_calculate_posterior_mean;
  int m_grid;
  bool m_disallow_empty_intervals_between_cps;
  bool m_write_cps_to_file;
  bool m_write_histograms_to_file;
  string m_model;
  bool m_importance_sampling;


 private:
  void usage(int status, char *);
  static struct option lopts[];

};



#endif
