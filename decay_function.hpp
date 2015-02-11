#ifndef DECAY_FN_H
#define DECAY_FN_H

#include "univariate_function.h"
#include <cmath>
#include <cfloat>

class Decay_Function : Univariate_Function
{
 public:
    Decay_Function(double rate = 1.0){ m_rate = rate; }
    virtual double function( double t ){ return exp(-m_rate*t); }
    virtual double cumulative_function( double t ){ return m_rate>0?(1-exp(-m_rate*t))/m_rate:t; }
    virtual double cumulative_function( double s, double t ) { return cumulative_function(t-s); }
 private:
    double m_rate;
};

#endif
