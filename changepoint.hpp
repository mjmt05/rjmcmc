#ifndef CHANGEPOINT_HPP
#define CHANGEPOINT_HPP

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
  double get_changepoint( changepoint* cpobj ){ return cpobj->getchangepoint(); } 

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





#endif
