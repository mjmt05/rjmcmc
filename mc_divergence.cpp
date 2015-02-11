#include "mc_divergence.hpp"

double* mc_divergence::phi_lookup = NULL;
double* mc_divergence::phi1_lookup = NULL;
double* mc_divergence::phi2_lookup = NULL;
unsigned int mc_divergence::length_phi_lookup = 0;
double* mc_divergence::chi_lookup = NULL;
unsigned int mc_divergence::length_chi_lookup = 0;
double* mc_divergence::log_lookup = NULL;
unsigned int mc_divergence::length_log_lookup = 0;
double mc_divergence::m_delta = 0.05;
double mc_divergence::m_phi_1=-gsl_sf_psi_int(1);


mc_divergence::mc_divergence(Divergence_Type divergence_type, Loss_Function loss_fn, unsigned int iterations){
  create_lookup_arrays(iterations,divergence_type,loss_fn);
  m_divergence_type = divergence_type;
  m_loss_fn = loss_fn;
  m_new_bin_divergence_change = divergence_type == BIAS ? m_phi_1 : 0;
  m_initial_sample_size = 1;
  m_initial_non_empty_bins = 1;
  reset();
}

void mc_divergence::reset(){
  m_divergence = 0;
  m_existing_bins_divergence = 0;
  m_new_bin_prob = 0;
  m_sample_size = 0;
  m_non_empty_bins = 0;
  m_waiting_time = m_last_waiting_time = 0;
}

mc_divergence::~mc_divergence(){
}

//  x/n -> (1-p)x/(n+1)+p y/(n+1)
//(n+1)x - (1-p)nx-pny
//  x+npx-npy
void mc_divergence::update_divergence( unsigned int n ){
  m_sample_size++;
  m_waiting_time++;
  if(n==1){//then non_empty_bins will have changed
    m_non_empty_bins++;
    m_last_waiting_time = m_waiting_time;
    m_waiting_time = 0;
  }
  if( m_loss_fn == AVERAGE && m_sample_size>m_initial_sample_size ){
    //      m_new_bin_prob = 1.0/(double)m_last_waiting_time;
    //      if( m_sample_size > m_last_waiting_time * m_non_empty_bins )
     m_new_bin_prob = (m_non_empty_bins-m_initial_non_empty_bins)/(double)(m_sample_size-m_initial_sample_size);
  }

  switch( m_divergence_type ){
  case BIAS:
    //    m_divergence += phi(n) - phi(n-1);
    if( m_loss_fn == MINIMAX )
      m_divergence += phi1(n-1);
    else{
      //      m_divergence -= n*phi2(n-1);
      m_existing_bins_divergence -= n*phi2(n-1);
      m_divergence = m_existing_bins_divergence;//not predict change in #bins
//      m_divergence = (1-m_new_bin_prob)*m_existing_bins_divergence + m_sample_size * m_new_bin_prob * m_new_bin_divergence_change;//this is wrong.
    }
    break;
  case CHI_SQUARE:
    if(n==1){
      m_existing_bins_divergence += m_new_bin_divergence_change;
	//      if( m_loss_fn == AVERAGE )
      m_new_bin_divergence_change = calculate_chi_squared_criteria(m_non_empty_bins+1) - m_existing_bins_divergence;
      if( m_loss_fn == MINIMAX )
	m_divergence = m_existing_bins_divergence;
    }
    if( m_loss_fn == AVERAGE )
//      m_divergence = m_existing_bins_divergence * m_sample_size;//not predict change in #bins
      m_divergence = m_existing_bins_divergence + m_sample_size * m_new_bin_prob * m_new_bin_divergence_change;
    break;
  case FOI:
    m_divergence = m_var;
    break;
  case EXTENT:
    m_divergence -= log1(n-1);
    break;
  case ENTROPY:
    m_divergence -= log1(n-1);
    break;
  case NONE:
    m_divergence = -n;
    break;
  }
}

double mc_divergence::get_divergence(){
  double divergence = m_divergence/m_sample_size;
  if( m_divergence_type == EXTENT || m_divergence_type == ENTROPY )
    divergence += log(m_sample_size);
  if( m_divergence_type == EXTENT)
    divergence = exp(2*divergence)/m_sample_size;

  if( m_loss_fn == MINIMAX )
    return divergence;
  return divergence/(m_sample_size+1);
}

double mc_divergence::phi( unsigned int n ){
    if( mc_divergence::phi_lookup && n < length_phi_lookup )
      return phi_lookup[ n ];
    else
      return (n>0?n*(log(n)-gsl_sf_psi_int(n)):0);// - pow(-1.0,(double)n)/(n+1.0);
}

double mc_divergence::phi1( unsigned int n ){
    if( mc_divergence::phi1_lookup && n < length_phi_lookup - 1 )
        return phi1_lookup[ n ];
    else
        return phi(n+1)-phi(n);
}

double mc_divergence::phi2( unsigned int n ){
    if( mc_divergence::phi2_lookup && n < length_phi_lookup - 2 )
        return phi2_lookup[ n ];
    else
        return phi(n+2)-2*phi(n+1)+phi(n);
}

void mc_divergence::create_phi_lookup( unsigned int max_n, Loss_Function loss_fn ){
  if(phi_lookup){
    if(length_phi_lookup >= max_n + 1)
      return;
    else
      delete [] phi_lookup;
  }
  unsigned int n;
  length_phi_lookup = max_n + 1;
  double* temp = new double[ length_phi_lookup ];
  temp[ 0 ] = 0;
  bool use_G = !true;
  if( use_G ){
    temp[ 1 ] = -log(2)+gsl_sf_psi_int(1);
    temp[ 2 ] = temp[ 1 ] + 2;
    for( n = 3; n < length_phi_lookup; n++ )
      if( (n%2) == 0 ) temp[ n ] = temp[ n - 2 ] + 2.0/(n-1);
      else temp[ n ] = temp[ n - 1 ];
    for( n = 1; n < length_phi_lookup; n++ )
      temp[n] = n*(log(n)-temp[n]);
  }
  else{
    for( n = 1; n < length_phi_lookup; n++ )
      temp[n] = phi(n);
  }
  phi_lookup = temp;
  if( loss_fn == MINIMAX ){
    phi1_lookup = new double[ length_phi_lookup - 1];
    for( n = 0; n < (length_phi_lookup-1); n++ )
      phi1_lookup[n]=phi_lookup[n+1]-phi_lookup[n];
  }
  if( loss_fn == AVERAGE ){
    phi2_lookup = new double[ length_phi_lookup - 2];
    for( n = 0; n < (length_phi_lookup-2); n++ )
      phi2_lookup[n]=phi_lookup[n+2]-2*phi_lookup[n+1]+phi_lookup[n];
    //	phi2_lookup[n]=phi1_lookup[n+1]-phi1_lookup[n];
    //    phi2_lookup[0]=phi(2)-2*phi(1);
  }
}

void mc_divergence::delete_phi_lookup(){
  if( phi_lookup )
    delete [] phi_lookup;
  if( phi1_lookup )
    delete [] phi1_lookup;
  if( phi2_lookup )
    delete [] phi2_lookup;
}

double mc_divergence::calculate_chi_squared_criteria(unsigned int n ){
  if( chi_lookup && n < length_chi_lookup )
    return chi_lookup[ n ];
  else
    return n>1?.5*gsl_cdf_chisq_Pinv(1-m_delta,n-1):0;
}

void mc_divergence::create_chi_lookup( unsigned int max_n ){
  if(chi_lookup){
    if(length_chi_lookup >= max_n + 1)
      return;
    else
      delete [] chi_lookup;
  }
  unsigned int n;
  length_chi_lookup = max_n + 1;
  double* temp = new double[ length_chi_lookup ];
  temp[ 0 ] = temp[ 1 ] = 0;
  for( n = 2; n <= max_n; n++ )
    temp[n] = calculate_chi_squared_criteria(n);
  chi_lookup = temp;
}

void mc_divergence::delete_chi_lookup(){
  if( chi_lookup )
    delete [] chi_lookup;
}

void mc_divergence::create_log_lookup( unsigned int max_n ){
  if(log_lookup){
    if(length_log_lookup >= max_n + 1)
      return;
    else
      delete [] log_lookup;
  }
  unsigned int n;
  length_log_lookup = max_n + 1;
  double* temp = new double[ length_log_lookup ];
  temp[ 0 ] = temp[ 1 ] = 0;
  for( n = 2; n <= max_n; n++ )
    temp[n] = n*log(n);
  log_lookup = temp;
}

void mc_divergence::delete_log_lookup(){
    if( log_lookup )
        delete [] log_lookup;
}

double mc_divergence::log1( unsigned int n ){
  if(!n)
    return 0;
  if( mc_divergence::log_lookup && n < length_log_lookup - 1 )
    return log_lookup[ n+1 ]-log_lookup[ n ];
  else
    return (n+1)*log(n+1)-n*log(n);
}


void mc_divergence::create_lookup_arrays(unsigned int length, Divergence_Type divergence_type, Loss_Function loss_fn)
{
  if(!length)
    return;
  if(divergence_type==BIAS)
    create_phi_lookup(length,loss_fn);
  else if(divergence_type==CHI_SQUARE)
    create_chi_lookup(length);
  else if(divergence_type==ENTROPY||divergence_type==EXTENT)
    create_log_lookup(length);
}

void mc_divergence::delete_lookup_arrays(Divergence_Type divergence_type)
{
  if(divergence_type==BIAS)
    delete_phi_lookup();
  else if(divergence_type==CHI_SQUARE)
    delete_chi_lookup();
  else if(divergence_type==ENTROPY||divergence_type==EXTENT)
    delete_log_lookup();
}

void mc_divergence::set_delta(double delta){
  m_delta=delta;
}

void mc_divergence::set_initial_number_of_bins(){
    m_initial_sample_size = m_sample_size;
    m_initial_non_empty_bins = m_non_empty_bins;
}
