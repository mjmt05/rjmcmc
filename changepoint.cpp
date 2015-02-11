#include "changepoint.hpp"

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

//
