#ifndef CHANGEPOINT_H
#define CHANGEPOINT_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <map>


 
using namespace std;


class changepoint{
  friend bool operator<(const changepoint &, const changepoint &);
  friend bool operator>(const changepoint &, const changepoint &);
  friend ostream &operator<< ( ostream& , const changepoint &);

 public:
  changepoint(double = 0,  int = 0, double = 0, double = 0);
  changepoint (const changepoint * f);
  ~changepoint();

  void setchangepoint(double xx){m_changepoint=xx;}
  void setdataindex(unsigned long long int yy){m_data_index=yy;}
  void setlikelihood(double zz){m_likelihood=zz;}
  void setmeanvalue(double mm){m_mean_value=mm;}
  void setvarvalue(double vv){m_var_value=vv;}
  void setdouble(double d){m_double=d;}
  double getchangepoint() const {return m_changepoint;}
  unsigned long long int getdataindex() const {return m_data_index;}
  double getlikelihood() const {return m_likelihood;}
  double getmeanvalue() const {return m_mean_value;}
  double getvarvalue() const {return m_var_value;}
  double getdouble() const{return m_double;}
  bool operator==(const changepoint *) const;
  void set_index_from_int(int value){m_data_index=m_vector_int[value];}
  int * get_int_vector() const {return m_vector_int;}
  double * get_double_vector() const {return m_vector_double;}
  int get_size_of_int() const{return m_size_of_int;}
  int get_size_of_double() const{return m_size_of_double;}
  void set_int_vector(int * vec, int size){m_vector_int = vec; m_size_of_int=size;} 
  void set_double_vector(double * dvec,int size){m_vector_double=dvec; m_size_of_double=size;}
  void set_general_pointer(map<unsigned int, unsigned int>* general_pointer){m_general_pointer = general_pointer;}
  map<unsigned int,unsigned int>* get_general_pointer() const { return m_general_pointer;}
  void delete_int_vector(){delete [] m_vector_int; m_size_of_int=0;}
  void delete_double_vector(){delete [] m_vector_double; m_size_of_double=0;}
   

 private:

  double m_likelihood,m_mean_value,m_var_value;
  double m_double;
  double * m_vector_double;
  int m_size_of_double;
  map<unsigned int, unsigned int> * m_general_pointer;
  
protected:
  double m_changepoint;
  unsigned long long int m_data_index;
  int * m_vector_int;
  int m_size_of_int;
};

changepoint::changepoint(double xx, int yy, double zz, double mm){
  setchangepoint(xx); setdataindex(yy); setlikelihood(zz); setmeanvalue(mm);
  m_size_of_int=0;
  m_var_value=0;
  m_size_of_double=0;
  m_general_pointer=NULL;
  m_vector_int = NULL;
  m_vector_double = NULL;
  m_double=0;

}

changepoint::changepoint(const changepoint * cp){
  m_changepoint = cp->m_changepoint;
  m_likelihood = cp->m_likelihood;
  m_mean_value = cp->m_mean_value;
  m_var_value = cp->m_var_value;
  m_data_index = cp->m_data_index;
  m_double = cp->m_double;
  m_size_of_int=cp->m_size_of_int;
  if(m_size_of_int>0){
    m_vector_int = new int[m_size_of_int];
    for(int i=0; i<m_size_of_int; i++)
      m_vector_int[i]=cp->m_vector_int[i];
    }

  m_size_of_double=cp->m_size_of_double;
  if(m_size_of_double>0){
    m_vector_double = new double[m_size_of_double];
    for(int i=0; i<m_size_of_double; i++)
      m_vector_double[i]=cp->m_vector_double[i];
  }

  if(cp->m_general_pointer)
    *m_general_pointer=*(cp->m_general_pointer);
  else 
    m_general_pointer=NULL;

}


changepoint::~changepoint(){
  if (m_general_pointer)
    delete m_general_pointer;

  if(m_size_of_int>0){
    delete [] m_vector_int;
  }
  if(m_size_of_double>0)
    delete [] m_vector_double;
}

bool operator<(const changepoint & cp1, const changepoint & cp2){
  return (cp1.getchangepoint()<cp2.getchangepoint());
}

bool operator>(const changepoint & cp1, const changepoint & cp2){
  return (cp1.getchangepoint()>cp2.getchangepoint());
}

ostream &operator<< ( ostream &output, const changepoint &cp){

    output << cp.m_changepoint<<" "<<cp.m_mean_value;
    return output;
}

bool changepoint::operator==(const changepoint * right) const{
  if(m_changepoint != right->m_changepoint)
    return false;
  if(m_mean_value != right->m_mean_value)
    return false;
  if(m_var_value !=right->m_var_value)
    return false;
  if(m_double != right->m_double)
    return false;
  if(m_size_of_int!=right->m_size_of_int)
    return false;
  if(m_size_of_int>0){
    for(int i=0; i<m_size_of_int; i++){
      if(m_vector_int[i]!=right->m_vector_int[i])
	return false;
    }
  }
  if(m_size_of_double!=right->m_size_of_double)
    return false;
  if(m_size_of_double>0){
    for(int i=0; i<m_size_of_double; i++){
      if(m_vector_double[i]!=right->m_vector_double[i])
	return false;
    }
  }
  return true;
 
}

double get_changepoint( changepoint* cpobj ){ return cpobj->getchangepoint(); }

#endif
