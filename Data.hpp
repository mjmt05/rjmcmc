#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <climits>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define FILENAMELENGTH 100

using namespace std;



template< class T>
class Data{


 public:
  Data(char *,unsigned long long int, unsigned long long int);
  Data(char *,bool,unsigned long long int, unsigned long long int);
  Data(string,bool=true);
  Data(istream& data_stream, unsigned int=0 );
  Data(T * , unsigned long long int=1, unsigned long long int=1, unsigned int =0);
 ~Data();
 
  Data<T>* get_subset(vector<unsigned int>* v);
  void read_from_file( bool = false, long long int = -1 );
  void from_store_to_data_matrix();
  void print_data_matrix();
  T** m_X;
  T** m_range_X;//the domain in which m_X lie.
  set<T> m_unique_data;
  unsigned long long int get_rows() const{return m_n;}
  unsigned long long int get_cols() const{return m_p;}
  unsigned long long int get_eventual_cols() const{return m_p_max;}
  T** get_range() const{return m_range_X;}
  T get_range_width()const{return m_range_X[0][1]-m_range_X[0][0];}
  T* copy_column(unsigned long long int);
  void data_construct();
  void construct(T * data, unsigned int offset = 0);
  bool is_empty() const{ return m_empty;}
  void increment_data_stream();
  T get_element( unsigned long long int, unsigned long long int );
  static void set_use_data_store( bool uds = true ){ m_use_data_store = uds; }
  unsigned long long int find_data_index(T , unsigned long long int = 0, unsigned long long int = 0, unsigned long long int = 0) const;
  void sort( unsigned int* ordering = NULL, unsigned int row = 0 );
  void order( unsigned int row = 0 );
  unsigned int* get_order(){ return m_order; }
  void kill_order(){ if(m_order) delete [] m_order; m_order = NULL; }
  void replace_with_modulo( T );
  void undo_replace_with_modulo();
  unsigned long long int find_unique_values();
  bool in_unqiue_values( T t ){ return m_unique_data.find(t) != m_unique_data.end(); }
  T*& operator[]( unsigned int row ){ return m_X[row]; }
  void replace_with_cumulative();
  void undo_replace_with_cumulative();
  void replace_with_right_cumulative();

 private:
  string m_filename;
  long long unsigned int m_n; //rows
  long long unsigned int m_p; //columns
  long long unsigned int m_p_max; //columns
  long long unsigned m_range_rows;
  bool m_empty;//set to true if the source data file is empty
  bool m_streaming;
  unsigned long long int m_window_size;//the number of columns held in the data matrix when m_streaming==true.
  unsigned long long int m_window_slide;//the number of columns incrementally read into the data matrix when m_streaming==true.
  vector<streamoff> m_stream_positions;//vector holding the position we have reached in each line of the data file when m_streaming==true.
  static bool m_use_data_store;
  unsigned int* m_order;
  vector<T> m_data_store;
  vector<unsigned int> m_modulo_store;
  T m_season;
};
template <class T>
bool Data<T>::m_use_data_store = false;

template <class T>
Data<T>::Data(char * file,unsigned long long int rows,unsigned long long int columns)
:m_filename(file), m_n(rows), m_p(columns)
{
  m_streaming = false;
  m_order = NULL;
  read_from_file(false,rows*columns);
  data_construct();
}

template <class T>
Data<T>::Data(string file,bool determine_rows)
:m_filename(file)
{
  m_n = 1;
  m_streaming = false;
  m_order = NULL;
  read_from_file( determine_rows );
  data_construct();
}

template <class T>
Data<T>::Data(char * file,bool determine_rows,unsigned long long int window_size, unsigned long long int window_slide )
:m_filename(file),m_window_size(window_size),m_window_slide(window_slide)
{
  m_n = 1;
  m_streaming = true;
  m_order = NULL;
  read_from_file( determine_rows );
  data_construct();
}


template <class T>
Data<T>::Data(istream& data_stream, unsigned int number_of_data ){
  T dummy;
  unsigned int counter = 0;
  while(data_stream >> dummy && (!number_of_data||++counter<number_of_data) ){
    m_data_store.push_back(dummy);
  }
  m_n = 1;
  m_p_max = m_p = m_data_store.size();
  T* data_array = new T[m_p];
  for(unsigned int i = 0; i < m_p; i++)
    data_array[i] = m_data_store[i];
  construct(data_array);
  delete [] data_array;
  m_data_store.clear();
}


template <class T>
Data<T>::Data(T * data, unsigned long long int rows, unsigned long long int columns, unsigned int offset)
:m_n(rows), m_p(columns), m_p_max(columns)
{
  construct(data,offset);
}

template <class T>
void Data<T>::construct(T * data, unsigned int offset)
{
  m_empty = true;
  m_streaming = false;
  m_order = NULL;
  unsigned int i;
//  m_filename=NULL;
  m_range_X=NULL;
  m_range_rows = m_p>1 ? m_n:1;
  if(m_n){
    m_X=new T* [m_n];
    m_range_X=new T* [m_range_rows];
    m_range_X[0]=new T[m_range_rows*2];
    if(m_p){
      m_empty = false;
      m_X[0]=new T[m_p*m_n];
      for (i=1;i<m_n;i++){
	m_X[i]=m_X[i-1]+m_p;
      }
      long long unsigned int row_counter=0;
      long long unsigned int col_counter=offset;
  
      for(i=0; i<m_p*m_n; i++){
	m_X[row_counter][col_counter]=data[i];

	unsigned long long int range_index = m_range_rows > 1 ? row_counter : 0;
	if( col_counter == 0 && range_index == row_counter)
	  m_range_X[range_index][0] = m_range_X[range_index][1] = m_X[row_counter][col_counter];
	else{
	  if(m_X[row_counter][col_counter]<m_range_X[range_index][0])
	    m_range_X[range_index][0] = m_X[row_counter][col_counter]; 
	  else if(m_X[row_counter][col_counter]>m_range_X[range_index][1])
	    m_range_X[range_index][1] = m_X[row_counter][col_counter];
	}
	if(col_counter==(m_p-1)) row_counter++;
	col_counter==(m_p-1)? col_counter=0 : col_counter++;
      }
    }
    else{ m_X[0]=new T[1];
      m_p=1;
    }
  }
  else{
    m_X=new T* [1];
    m_X[0]=new T[1];
  }
}

template <class T>
Data<T>* Data<T>::get_subset(vector<unsigned int>* v){
  unsigned int v_size = v->size();
  T* data = new T[ m_n * v_size ];
  for(unsigned int i = 0; i< m_n; i++ )
    for(unsigned int j = 0; j< v_size; j++ )
      data[ i * m_n + j ] = m_X[ i ][ (*v)[j] ];
  Data<T>* data_obj = new Data( data, m_n, v_size );
  delete [] data;
  return data_obj;
}

template <class T>
Data<T>::~Data()
{
  if(m_X){
    delete[] m_X[0];
    delete [] m_X;
    m_X=0;
  }
  if(m_range_X){
    delete[] m_range_X[0];
    delete [] m_range_X;
    m_range_X=0;
  }
  

}

template <class T>
void Data<T>::data_construct(){
  unsigned int i;
  if(m_n){
    m_range_rows = m_p>1?m_n:1;
    m_X=new T* [m_n];
    m_range_X=new T* [m_range_rows];
    m_range_X[0]=new T[m_range_rows*2];
    if(m_p) m_X[0]=new T[m_p*m_n];
    else{ m_X[0]=new T[1];
      m_p=1;
      m_X[0][0]=0;//1e7;
      m_range_X[0][0] = m_range_X[0][1] = m_X[0][0];
    }
  }
  else{
    m_range_rows = 1;
    m_X=new T* [1];
    m_X[0]=new T[1];
    m_X[0][0]=0;//1e7;
    m_range_X=new T* [m_range_rows];
    m_range_X[0]=new T[m_range_rows*2];
    m_range_X[0][0] = m_range_X[0][1] = m_X[0][0];
  }

  for (i=1;i<m_n;i++){
    m_X[i]=m_X[i-1]+m_p;
    if(m_range_rows>1)
      m_range_X[i]=m_range_X[i-1]+2;
  }
 
  from_store_to_data_matrix();
}


template <class T>
void Data<T>::read_from_file( bool determine_number_of_rows, long long int how_many )
{
  T filecontents;
  long long int total_count = 0, row_count = 0;
  string line;
  ifstream inFile(m_filename.c_str());
  if (!inFile.is_open()){
    m_empty = true;
    m_n = m_p = 0;
    m_X = NULL;
    m_range_X = NULL;
    //    return;
    cerr << "Error: " <<  m_filename << " could not be opened." <<endl;
    exit(1);
  }

  m_p = 0;
  while( inFile.good() ){
    getline(inFile,line);
    if(line.size()>0){
      unsigned int col_count = 0;
      row_count++;
      istringstream iss(line); 
      while( iss >> filecontents && ( how_many < 0 || total_count++ < how_many)  ){
	col_count++;
	if(m_use_data_store && ( !m_streaming || col_count <= m_window_size ) )
	  m_data_store.push_back( filecontents );
	if(m_streaming && (col_count == m_window_size) ){
	  streamoff offset = iss.tellg();
	  m_stream_positions.push_back(offset);
      }
	m_p++;
      }
    }
  }
  if( determine_number_of_rows )
    m_n = row_count;
  if((long long int )m_p < how_many ){
    cerr << "Error: The file " << m_filename << " does not contain " << how_many << " values." << endl;
    exit(1);
  }
  if(!m_p){
    m_empty = true;
    m_n = 0;
  }else{
    m_empty = false;
    m_p /= m_n;
  }
  m_p_max = m_p;
  if(m_streaming)
    m_p = m_window_size;
  inFile.close();
}

template <class T>
void Data<T>::from_store_to_data_matrix()
{
  unsigned int counter = 0;
  ifstream inFile;
  if(!m_use_data_store)
      inFile.open(m_filename.c_str());
  bool next_line_with_next_row = m_p_max > m_p;
  for (unsigned int i=0; i<m_n; i++){
    string line;
    istringstream iss;
    if(!m_use_data_store&&next_line_with_next_row){
      getline(inFile,line);
      iss.str(line); 
    }
    for (unsigned int j=0; j<m_p; j++){
      if(m_use_data_store)
	m_X[i][j] = m_data_store[counter++];
      else{
	if(next_line_with_next_row)
	  iss >> m_X[i][j];
	else
	  inFile >> m_X[i][j];
      }
      if(m_range_X){
	unsigned int range_index = m_range_rows > 1 ? i : 0;
	if( j == 0 && range_index == i)
	  m_range_X[range_index][0] = m_range_X[range_index][1] = m_X[i][j];
	else{
	  if(m_X[i][j]<m_range_X[range_index][0])
	    m_range_X[range_index][0] = m_X[i][j]; 
	  else if(m_X[i][j]>m_range_X[range_index][1])
	    m_range_X[range_index][1] = m_X[i][j];
	}
      }
    }
  }
  if(m_use_data_store)
    m_data_store.clear();
  else
    inFile.close();
}

template <class T>
void Data<T>::increment_data_stream(){
  unsigned int i=0, j=0;
  for(i=0; i<m_n; i++)
    for(j=0; j<(m_window_size-m_window_slide); j++)
      m_X[i][j] = m_X[i][j+m_window_slide];
  m_p += m_window_slide;

  ifstream inFile(m_filename.c_str());
  string line;
  i=0;
  T filecontents;
  while( inFile.good() ){
    getline(inFile,line);
    if(line.size()>0){
      j = 0;
      istringstream iss(line);
      iss.seekg(m_stream_positions[i]);
      while( j < m_window_slide && iss >> filecontents ){
	  m_X[i][m_window_size-m_window_slide+j] = filecontents;
	  j++;
      }
      m_stream_positions[i]=iss.tellg();
      i++;
    }
  }
}

template <class T>
T Data<T>::get_element( unsigned long long int i , unsigned long long int j ){
  if(!m_streaming)
    return m_X[i][j];
  return m_X[i][j-(m_p-m_window_size)];
}

template <class T>
T* Data<T>::copy_column(unsigned long long int col ){
  if(col>=m_p)
    return NULL;
  T* v = new T[ m_n ];
  for(unsigned int i=0; i<m_n; i++)
    v[i] = m_X[ i ][ col ];
  return v;
}

template <class T>
void Data<T>::print_data_matrix()
{
  for (unsigned int i=0; i<m_n; i++){
    for (unsigned int j=0; j<(m_streaming?m_window_size:m_p); j++)
      cout << m_X[i][j] << " ";
    cout << endl;
  }
}


template <class T>
unsigned long long int Data<T>::find_data_index(T theta, unsigned long long int i, unsigned long long int l, unsigned long long int h) const
{
  if(m_empty)
    return 0;//static_cast<unsigned int>(theta);
  if(l == m_p)
    return m_p;
  if(is_empty()||m_X[i][l]>=theta)
    return l;
  if(!h)
    h = m_p - 1;
  if(m_X[i][h]<theta)
    return h+1;
  
  unsigned long long int m = (l+h)/2;
  while( h-l > 1 ){
    ( m_X[i][m]>=theta ) ? h = m : l = m;
    m = (l+h)/2;
  }
  return m+1;
}

template <class T>
void Data<T>::order( unsigned int row ){
  if(m_empty)
    return;
  pair<T,unsigned int>* v = new pair<T,unsigned int>[m_p];
  for( unsigned int j = 0; j < m_p; j++ )
    v[j] = make_pair(m_X[row][j],j);
  std::sort(v,v+m_p);
  if( m_order )
    delete [] m_order;
  m_order = new unsigned int[m_p];
  for( unsigned int j = 0; j < m_p; j++ )
    m_order[j] = v[j].second;
  //    cout << m_order[j] << endl;
  delete [] v;
}


template <class T>
void Data<T>::sort( unsigned int* ordering, unsigned int row ){
  if(m_empty)
    return;
  if(!ordering){
    order(row);
    ordering = m_order;
  }
  T** old_mx = m_X;
  m_X=new T*[m_n];
  m_X[0]=new T[m_p*m_n];
  unsigned int i,j;
  for(i=1;i<m_n;i++)
    m_X[i]=m_X[i-1]+m_p;
  for(i=0;i<m_n;i++)
    for(j=0;j<m_p;j++)  
      m_X[i][j] = old_mx[i][ordering[j]];
  delete [] old_mx[0];
  delete [] old_mx;  
}

template <class T>
unsigned long long int Data<T>::find_unique_values(){
  m_unique_data.clear();
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned int j = 0; j < m_p; j++ )
      m_unique_data.insert(m_X[i][j]);
  return m_unique_data.size();
}

template <class T>
void Data<T>::replace_with_modulo( T val ){
  m_season = val;
  unsigned int counter = 0;
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned int j = 0; j < m_p; j++ ){
      m_modulo_store.push_back(static_cast<unsigned int>(m_X[i][j]/val));
      m_X[i][j] -= m_modulo_store[counter++] * val;
    }
}

template <class T>
void Data<T>::undo_replace_with_modulo(){
  unsigned int counter = 0;
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned int j = 0; j < m_p; j++ )
      m_X[i][j] += m_modulo_store[m_order?m_order[counter++]:counter++] * m_season;
}

template <class T>
void Data<T>::replace_with_cumulative(){
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned int j = 1; j < m_p; j++ )
      m_X[i][j] += m_X[i][j-1];
}

template <class T>
void Data<T>::undo_replace_with_cumulative(){
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned long long int j = m_p-1; j >0; j-- )
      m_X[i][j] -= m_X[i][j-1];
}

template <class T>
void Data<T>::replace_with_right_cumulative(){
  for( unsigned int i = 0; i < m_n; i++ )
    for( unsigned int j = 1; j < m_p; j++ )
      m_X[i][m_p-j-1] += m_X[i][m_p-j];
}

#endif
