#ifndef UV_FN_HPP
#define UV_FN_HPP

class Univariate_Function
{
 public:
  virtual double function( double t ) = 0;
  virtual double cumulative_function( double t ) = 0;
  virtual double cumulative_function( double s, double t ) = 0;
  virtual ~Univariate_Function(){;}

};


#endif

