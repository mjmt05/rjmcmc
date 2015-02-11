#ifndef UV_FN_H
#define UV_FN_H

class Univariate_Function
{
 public:
  virtual double function( double t ) = 0;
  virtual double cumulative_function( double t ) = 0;
  virtual double cumulative_function( double s, double t ) = 0;
  virtual ~Univariate_Function(){;}

};


#endif

