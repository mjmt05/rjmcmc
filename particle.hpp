#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>

/*class assumes that the array you pass in to create an object is already sorted*/

using namespace std;

template<class T> class Particle;

template<class T>
ostream& operator<<(ostream&,const Particle<T> &);

template< class T>
class Particle
{
  friend ostream &operator<< <> ( ostream& , const Particle<T> &);

 public:
  Particle(int=0 ,T ** = NULL, T* = NULL, unsigned int birth_time=0 ); 
  Particle(const Particle *, const Particle *, unsigned int birth_time=0);
  ~Particle();
  void settheta(T**);
  void delete_component(unsigned int);
  void add_component(T*,unsigned int);
  void change_component(T*,int);
  unsigned int get_dim_theta() const {return m_dim_theta;} 
  T* get_theta_component(int k) const;
  T* get_last_theta_component() const {return get_theta_component(m_dim_theta-1);}
  unsigned int find_position(T *, bool = false, unsigned int = 0, unsigned int = 0);
  void set_log_posterior(double post) {m_log_posterior = post;}
  double get_log_posterior() const {return m_log_posterior;}
  T* get_intercept() const{return m_intercept;}
  bool operator==(const Particle<T> *) const;
  void set_weight( double weight ){ m_log_weight = weight; }
  long double get_weight(){ return m_log_weight; }
  unsigned int get_birth_time(){ return m_birth_time; }
  bool does_particle_exist(T *);


 protected:

  unsigned int m_dim_theta;
  T **m_theta;
  T * m_intercept;
  double m_log_posterior;
  long double m_log_weight;
  unsigned int m_birth_time;
  static int particleCount;
  void sort( T **, unsigned int);
  void swap( T * const, T * const);

};

template<class T> 
int Particle<T>::particleCount=0;


template <class T>
Particle<T>::Particle(int k,T ** thetaarray, T* thetaintercept, unsigned int birth_time)
:m_dim_theta(k)
{
  settheta(thetaarray);

  m_intercept = thetaintercept;
  
  m_birth_time=birth_time;
  ++ particleCount;
}


template <class T>
Particle<T>::Particle(const Particle *particle1, const Particle *particle2, unsigned int birth_time)
{

  if(particle2==NULL){
    m_dim_theta=particle1->m_dim_theta;
    m_log_weight = particle1->m_log_weight;

    if (particle1->m_intercept != NULL){
      m_intercept = new T(particle1->m_intercept);
    }
    else{
      m_intercept = NULL;
    }

    if (m_dim_theta==0){
      m_theta=NULL;
    }
    else{
      m_theta=new T*[m_dim_theta];
      assert(m_theta != 0);
                    
      for (unsigned int i=0; i<particle1->m_dim_theta; i++){
	m_theta[i] = new T(particle1->m_theta[i]);
      }

    }
    m_log_posterior = particle1->m_log_posterior;
  }else{
    if (particle1->m_intercept!=NULL && particle2->m_intercept !=NULL){
      m_intercept = new T;
            
      if (*(particle1->m_intercept)< *(particle2->m_intercept)){
	*m_intercept = *(particle1->m_intercept);
      }
       
      else{
	*m_intercept = *(particle2->m_intercept);
      }
    }
    else{
      m_intercept = NULL;
    }
        
        
    int k1=particle1->m_dim_theta;
    int k2=particle2->m_dim_theta;
    m_dim_theta=k1+k2;
    m_log_weight = particle1->m_log_weight + particle2->m_log_weight;
	
    if (m_dim_theta >  0)
      {
	m_theta=new T*[m_dim_theta];
	assert(Particle<T>::m_theta != 0);

	for (int i=0; i<k1; i++)
	  {
	    m_theta[i] = new T;
	    *(m_theta[i])=*(particle1->m_theta[i]);
	  }

	for (int i=k1; i<(k2+k1); i++)
	  {
	    m_theta[i]= new T;
	    *(m_theta[i])=*(particle2->m_theta[i-k1]);
	  }
      }
    else{
      m_theta=NULL;
    }

    m_log_posterior = particle1->m_log_posterior + particle2->m_log_posterior;

  }
  
  if(particle1->m_birth_time>birth_time)
    m_birth_time=particle1->m_birth_time;
  else
    m_birth_time=birth_time;
  
  ++particleCount;
}

template <class T>
Particle<T>::~Particle()
{
  for (unsigned int i=0; i<m_dim_theta;i++){
    delete m_theta[i];
  }

  if (m_dim_theta>0) {
    delete [] m_theta;
  }

  if(m_intercept){
    delete m_intercept;
  }
      
  --particleCount;
}


template <class T>
void Particle<T>::settheta(T ** thetaarray)
{

  if (m_dim_theta==0){
    m_theta=NULL;
  }else{
    m_theta=thetaarray;
  }	

}

template <class T>
void Particle<T>::delete_component(unsigned int l)
{

  if (m_dim_theta==0){
    cerr << "In Particle.h: the dimension of the particle is 0, there is no objects to delete" << endl;
  }

  T ** temp_theta=m_theta;

  if (m_dim_theta-1==0){
    m_theta=NULL;
  }else{
    m_theta = new T*[m_dim_theta-1];
  }
     
  for (unsigned int i=0; i<l; i++){
    m_theta[i]=temp_theta[i];
  }
  
  for (unsigned int i=l; i<m_dim_theta-1; i++){
    m_theta[i]=temp_theta[i+1];
  }

  --m_dim_theta;

  delete temp_theta[l];
  delete [] temp_theta;

}


template <class T>
void Particle<T>::add_component(T* new_component, unsigned int l)
{
  
  T** temp_theta = m_theta;
  m_theta = new T*[m_dim_theta+1];

  for (unsigned int i=0; i<l; i++){
    m_theta[i]=temp_theta[i];
  }

  m_theta[l] = new_component;

  for(unsigned int i=l+1; i<m_dim_theta+1; i++){
    m_theta[i]=temp_theta[i-1];
  }

  if (m_dim_theta>0){
    delete[] temp_theta;
  }

  ++m_dim_theta;
}

template <class T>
void Particle<T>::change_component(T* new_value,int l)
{
  if(l<0){
    delete m_intercept;
    m_intercept=new_value;
  }else{
    delete m_theta[l];
    m_theta[l]=new_value;
  }

}

template<class T>
void Particle<T>::sort(T **arraytosort,unsigned int size)
{
  for (int i=0; i<size-1; i++)
    for (int j=0; j<size-1; j++)
      if( *arraytosort[j] > *arraytosort[j+1] )
	swap( arraytosort[j], arraytosort[j+1] );
}

template<class T>
void Particle<T>::swap( T * const element1Ptr, T * const element2Ptr )
{
  T hold = *element1Ptr;
  *element1Ptr = *element2Ptr;
  *element2Ptr = hold;
}


template <class T>
ostream &operator<< ( ostream &output, const Particle<T> &p)
{

  output << *p.m_intercept << ' ';
  for (unsigned int i=0; i<p.m_dim_theta; i++){
    output << *p.m_theta[i] << ' ';
  }
  output<<endl;
  return output;
}

template<class T>
unsigned int Particle<T>::find_position(T * new_theta, bool bisection,  unsigned int l, unsigned int h ){ 
  
  if( l == m_dim_theta )
    return m_dim_theta;
  if(m_dim_theta==0||*new_theta<*m_theta[l])
    return l;
  if(!h)
    h = m_dim_theta-1;
  if(*new_theta>*m_theta[h])
    return h+1;
  if( bisection ){
    unsigned int m = (l+h)/2;
    while( h-l > 1 ){
      ( *new_theta<*m_theta[m] ) ? h = m : l = m;
      m = (l+h)/2;
    }
    return m+1;
  }
  
  for(unsigned int j=l; j<m_dim_theta; j++)
    if(*new_theta<*m_theta[j])
      return j;
  

  /*a location for new_theta could not be found*/
  cerr<<"Particle.h: a location for the new object in the particle could not be found " << endl;
  cerr<<"New object = "<<*new_theta<<endl;
  for (unsigned int j =0; j < m_dim_theta; j++) {
    cerr<<"Particle = "<<*m_theta[j]<<" ";
  }
  cerr << endl;
  
  return m_dim_theta +1;
}

template<class T>
bool Particle<T>::does_particle_exist(T* new_theta) {
  for (unsigned int j = 0; j < m_dim_theta; j++) {
    if (new_theta->getchangepoint() == m_theta[j]->getchangepoint()) {
      return true;
    }
  }
  return false;

}
  
template<class T>
  T* Particle<T>::get_theta_component(int k)const {

  if (k<0)
    return m_intercept;
  else
    return m_theta[k];
}

template<class T>
bool Particle<T>::operator==(const Particle<T> *right) const{

  if(m_dim_theta!=right->m_dim_theta)
    return false;


  for(int i=0; i<m_dim_theta; i++)
    if (*m_theta[i]!=*(right->m_theta[i]))
      return false;


  return true;
}


#endif


