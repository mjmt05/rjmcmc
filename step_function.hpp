#ifndef STEPFN_HPP
#define STEPFN_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <vector>
#include "univariate_function.hpp"
#include "Data.hpp"
#define GSL_RANGE_CHECK_OFF

using namespace std;

class Step_Function : public Univariate_Function
{

public:
  Step_Function(double* knots, double* heights, unsigned int num_knots, double end = DBL_MAX, bool do_cumulative = true);
  Step_Function(vector<double>* knots, vector<double>* heights, double end = DBL_MAX, bool do_cumulative = true);
  Step_Function(const string input_filename, double end = DBL_MAX, bool do_cumulative = true);
  Step_Function( const Step_Function* f );
  ~Step_Function();
  void construct_step_function(vector<double>* knots, vector<double>* heights, double end, bool do_cumulative);
  void calculate_cumulative_function();
  unsigned int find_position_in_knots( double t, bool bisection = false, unsigned int l = 0, unsigned int h = 0 );
  virtual double function( double t );
  virtual double cumulative_function( double t );
  virtual double cumulative_function( double s, double t );
  double get_end_cumulative(){ return m_end_cumulative; }
  double* get_changepoints(){ return m_knots; }
  double* get_heights(){ return m_heights; }
  double* get_cumulative_heights(){ return m_cumulative; }
  unsigned int get_num_changepoints(){ return m_num_knots; }
  Step_Function* combine_with_step_function( Step_Function* f, bool multiply = true );
  Step_Function* add_step_function( Step_Function* f ){ return combine_with_step_function( f, false ); };
  Step_Function* multiply_step_function( Step_Function* f ){ return combine_with_step_function( f, true ); };
  void scale_step_function(double);
  void standardise_step_function(){ scale_step_function(1/m_end_cumulative); }
  double get_end(){ return m_end; }
  friend ostream& operator<< ( ostream& os, Step_Function* & f );
  vector<unsigned int>* find_data_segments( Data<double>* D );

protected:
  double* m_knots;//assume in increasing order
  double* m_heights;
  double m_end;//after this value, the step function gets repeated, infinitely.
  double m_end_cumulative;//the integral of the step function up to m_end
  double* m_cumulative;//the integral of the step function up to each knot point
  unsigned int m_num_knots;
};




#endif
