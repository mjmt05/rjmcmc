#ifndef STEPFN_H
#define STEPFN_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cassert>
#include <vector>
#include "univariate_function.h"
#include "Data.h"
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

Step_Function::Step_Function(double* knots, double* heights, unsigned int num_knots, double end, bool do_cumulative)
{
  m_num_knots = num_knots;
  m_end = end;
  m_knots = new double[ m_num_knots ];
  m_heights = new double[ m_num_knots ];
  for(unsigned int i = 0; i < m_num_knots; i++ ){
    m_knots[ i ] = knots[ i ];
    m_heights[ i ] = heights[ i ];
  }
  m_end_cumulative = 1; 
  m_cumulative = NULL;
  if(do_cumulative)
    calculate_cumulative_function();
}

Step_Function::Step_Function(vector<double>* knots, vector<double>* heights, double end, bool do_cumulative)
{
   construct_step_function(knots,heights,end,do_cumulative);
}

Step_Function::Step_Function(const string input_filename, double end, bool do_cumulative)
{
  m_num_knots = 0;
  ifstream InputStream(input_filename.c_str(), ios::out);
  if(!InputStream.is_open())
    return;
  vector<double> knots, heights;
  double knot, height;
  string line;
  while(!InputStream.eof() ){
    getline(InputStream,line);
    if(line.size()>0){
      istringstream iss(line); 
      if(!(iss >> knot >> height))
	break;
      knots.push_back(knot);
      heights.push_back(height);
      m_num_knots++;
    }
  }
  if(!m_num_knots)
    return;
  construct_step_function(&knots,&heights,end,do_cumulative);
}

Step_Function::Step_Function( const Step_Function* f )
{
  m_num_knots = f->m_num_knots;
  if(m_num_knots){
    m_knots = new double[ m_num_knots ];
    m_heights = new double[ m_num_knots ];
    for(unsigned int i = 0; i < m_num_knots; i++ ){
      m_knots[ i ] = f->m_knots[ i ];
      m_heights[ i ] = f->m_heights[ i ];
    }
  }
  m_cumulative = NULL;
  if(f->m_cumulative)
    calculate_cumulative_function();
  else
    m_end_cumulative = 1;
}

void Step_Function::construct_step_function(vector<double>* knots, vector<double>* heights, double end, bool do_cumulative){
  m_num_knots = knots->size();
  m_end = end;
  if( m_num_knots ){
    m_knots = new double[ m_num_knots ];
    m_heights = new double[ m_num_knots ];
    for(unsigned int i = 0; i < m_num_knots; i++ ){
      m_knots[ i ] = (*knots)[ i ];
      m_heights[ i ] = (*heights)[ i ];
    }
  }
  m_cumulative = NULL;
  if(do_cumulative)
    calculate_cumulative_function();
  else
    m_end_cumulative = 1; 
}

Step_Function::~Step_Function(){
  if( m_num_knots ){
    delete [] m_knots;
    delete [] m_heights;
    if(m_cumulative)
      delete [] m_cumulative;
  }
}

void Step_Function::calculate_cumulative_function(){
  m_cumulative = new double[m_num_knots];
  m_cumulative[0] = 0;
  for(unsigned int i = 1; i < m_num_knots; i++ ){
    m_cumulative[i] = m_cumulative[i-1]+m_heights[i-1]*(m_knots[i]-m_knots[i-1]);
  }
  if(m_end<DBL_MAX)
    m_end_cumulative = m_cumulative[m_num_knots-1]+m_heights[m_num_knots-1]*(m_end-m_knots[m_num_knots-1]);
}

unsigned int Step_Function::find_position_in_knots( double t, bool bisection, unsigned int l, unsigned int h ){
  if( l == m_num_knots )
    return m_num_knots;
  if(m_num_knots==0||t<m_knots[l])
      return l;
    if(!h)
      h = m_num_knots-1;
    if(t>m_knots[h])
      return h+1;
    if( bisection ){
      unsigned int m = (l+h)/2;
      while( h-l > 1 ){
	( t<m_knots[m] ) ? h = m : l = m;
	m = (l+h)/2;
      }
      return m+1;
    }

    for(unsigned int j=l; j<m_num_knots; j++)
      if(t<m_knots[j])
	return j;
    return m_num_knots;
}

double Step_Function::function( double t ){
  if(t>m_end){
    unsigned int num_repeats = (unsigned int)(t/m_end);
     t -= num_repeats * m_end;
  }
  return m_heights[ find_position_in_knots(t) - 1 ];
}

double Step_Function::cumulative_function( double t ){
  unsigned int num_repeats = 0;
  if(t>m_end){
    num_repeats = (unsigned int)(t/m_end);
     t -= num_repeats * m_end;
  }
  unsigned int t_pos = find_position_in_knots(t);
  return num_repeats * m_end_cumulative + m_cumulative[t_pos-1] + (t-m_knots[t_pos-1])*m_heights[t_pos-1];
}

double Step_Function::cumulative_function( double s, double t ){
  unsigned int num_repeats_s = 0, num_repeats_t = 0;
  if(s>m_end){
    num_repeats_s = (unsigned int)(s/m_end);
     s -= num_repeats_s * m_end;
  }
  if(t>m_end){
    num_repeats_t = (unsigned int)(t/m_end);
     t -= num_repeats_t * m_end;
  }
  double cum = (num_repeats_t - num_repeats_s) * m_end_cumulative;
  unsigned int s_pos = find_position_in_knots(s);
  unsigned int t_pos = find_position_in_knots(t);//,true,s_pos);//assumes t>s
  if(s_pos==t_pos)
    return cum+(t-s)*m_heights[t_pos-1];
  return cum+m_cumulative[t_pos-1] + (t-m_knots[t_pos-1])*m_heights[t_pos-1] - m_cumulative[s_pos-1] - (s-m_knots[s_pos-1])*m_heights[s_pos-1];
}

Step_Function* Step_Function::combine_with_step_function( Step_Function* f, bool multiply ){
  if(!f->m_num_knots)
    return new Step_Function(this);
  unsigned int counter1 = 0, counter2 = 1;//assume first knots the same
  double season_intercept2 = 0;
  vector<double> new_knots;
  vector<double> new_heights;
  double new_end = m_end;
  if(f->m_end>new_end)
    new_end = f->m_end;
  double f1 = 0;
  if(m_heights)
    f1 = m_heights[0];
  double f2 = f->m_heights[0];
  while( counter1<m_num_knots || counter2<f->m_num_knots ){
    if( counter1 == 0 || (counter2 >= f->m_num_knots) || ( (counter1 < m_num_knots) && (m_knots[ counter1 ] <= season_intercept2 + f->m_knots[ counter2 ]) ) ){
      new_knots.push_back(m_knots[counter1]);
      f1 = m_heights[counter1];
      counter1++;
    }
    else{
      new_knots.push_back(season_intercept2 + f->m_knots[counter2]);
      f2 = f->m_heights[counter2];
      counter2++;
      if(counter2 == f->m_num_knots && f->m_end<DBL_MAX){
	season_intercept2 += f->m_end;
	bool keep_going=false;
	if(m_end<DBL_MAX)
	  keep_going = season_intercept2+f->m_knots[0]<new_end;
	else
	  keep_going = counter1<m_num_knots;
	if(keep_going)
	  counter2=0;
      }
    }
    new_heights.push_back(multiply ? f1*f2 : f1+f2);
  }
  Step_Function* new_step_function = new Step_Function( &new_knots, &new_heights, new_end, m_cumulative!=NULL );
  return new_step_function;
}

void Step_Function::scale_step_function(double a){
  for(unsigned int i = 0; i < m_num_knots; i++ ){
    m_heights[ i ] *= a;
    if(m_cumulative)
      m_cumulative[ i ] *= a;
  }
  m_end_cumulative *= a;
}

ostream& operator<< ( ostream& os, Step_Function* & f ){
  if(!f){
    os << 0 << endl;
    return os;
  }
  os.setf( ios::fixed );
  for(unsigned int i = 0; i < f->m_num_knots; i++ ){
    os << f->m_knots[ i ] << "\t" << f->m_heights[ i ];
    if(f->m_cumulative)
      os << "\t" << f->m_cumulative[ i ];
    os<< endl;
  }
  return os;
}

vector<unsigned int>* Step_Function::find_data_segments( Data<double>* D ){
  unsigned int num_season_segments = m_num_knots + 1;
  vector<unsigned int>* m_data_seasons = new vector<unsigned int>[ num_season_segments ];
  unsigned int num_cols = D->get_cols();
  for(unsigned int i = 0; i < num_cols; i++){
    double t = (*D)[0][i];
    t = t - m_end * static_cast<unsigned int>(t/m_end);
    unsigned int j = 0;
    while(j<num_season_segments && t>m_knots[j] )
      j++;
    m_data_seasons[j].push_back(i);
  }
  return m_data_seasons;
}



#endif
