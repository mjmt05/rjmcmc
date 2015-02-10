#ifndef HISTOGRAM_TYPE_H
#define HISTOGRAM_TYPE_H

#include "histogram.h"
#include "particle.h"

template<class T>
double double_dereference(T* t){ return *static_cast<double*>(static_cast<void*>(t)); }

template< class T>
class Histogram_Type : public Histogram
{
  typedef double (*Func_Ptr)(T*);

 public:
  Histogram_Type(double start, double end, unsigned int num_bins, double m_bin_width=-1, bool bounded=true, bool weighted=false, bool one_d=false, unsigned int max_dim=0, string* input_filename=NULL, bool calculate_divergence=false, Divergence_Type divergence_tgype=BIAS, Loss_Function loss_fn=MINIMAX, unsigned int look_up_length=0,bool estimate_autocorrelation=true);
  Histogram_Type( const Histogram_Type<T>& h );
  ~Histogram_Type();
  void set_value_function( Func_Ptr get_component_value ){ m_get_component_value = get_component_value; }
  void calculate_bin( Particle<T>* particle, unsigned int=0 );
  void use_single_component(unsigned short int option = 1){ m_use_single_component = option; }
  
 protected:
  Func_Ptr m_get_component_value;
  unsigned short int m_use_single_component;
  vector<unsigned int> m_previous_bin;
  unsigned int m_previous_dim;//the dimension of the previous particle
};

template <class T>
Histogram_Type<T>::Histogram_Type(double start, double end, unsigned int num_bins, double bin_width, bool bounded, bool weighted, bool one_d,unsigned int max_dim, string* input_filename, bool calculate_divergence, Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int look_up_length,bool estimate_autocorrelation):Histogram(start,end,num_bins,bin_width,bounded,weighted,one_d,max_dim,input_filename,calculate_divergence,divergence_type,loss_fn,look_up_length,estimate_autocorrelation){
  m_get_component_value = &double_dereference;//NULL;
  m_use_single_component = 0;
}

template <class T>
Histogram_Type<T>::Histogram_Type(const Histogram_Type<T>& h ):Histogram(h){
  m_get_component_value = h.m_get_component_value;
  m_use_single_component = 0;
  m_estimate_autocorrelation = h.m_estimate_autocorrelation;
}


template <class T>
Histogram_Type<T>::~Histogram_Type(){}

template <class T>
void Histogram_Type<T>::calculate_bin( Particle<T>* particle, unsigned int particle_dim ){
  if(m_estimate_autocorrelation && m_samples > 0){
      m_previous_bin = m_current_bin;
      m_previous_dim = m_current_dim;
  }
  m_current_bin.clear();
  m_current_dim = 0;
  bool same_particle_as_last_time = m_estimate_autocorrelation && (m_samples > 0) && (m_previous_dim == particle_dim);
  if(!m_use_single_component)
      for(unsigned int k=0; k<particle_dim; k++){
          T* component_k = particle->get_theta_component(k);
          double val = m_get_component_value( component_k );//m_get_component_value? m_get_component_value( component_k ): static_cast<double>(*component_k);
          Histogram::calculate_bin( val );
          if(same_particle_as_last_time)
            same_particle_as_last_time = m_previous_bin[k] == m_current_dim_bin;
      }
  else if(particle_dim>0){
      if(m_use_single_component==2)
          Histogram::calculate_bin( m_get_component_value( particle->get_theta_component(particle_dim-1) ) );
      else if(m_use_single_component==1 || m_use_single_component==3)
          Histogram::calculate_bin( m_get_component_value( particle->get_theta_component(0) ) );
      if(m_use_single_component==3 && particle_dim>1)
          Histogram::calculate_bin( m_get_component_value( particle->get_theta_component(particle_dim-1) ) );
      if(m_use_single_component==4){
	double old_val = m_start;
	double max_dist = -1;
	for(unsigned int k=0; k<particle_dim; k++){
          T* component_k = particle->get_theta_component(k);
          double val = m_get_component_value( component_k );//m_get_component_value? m_get_component_value( co
	  if(val-old_val>max_dist)
	    max_dist = val-old_val;
	  old_val = val;
	}
	if(m_end-old_val>max_dist)
	  max_dist = m_end-old_val;
	Histogram::calculate_bin( max_dist );
      }
  }
  if(same_particle_as_last_time)
    m_num_bin_repeats++;
}


#endif
