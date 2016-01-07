#include "Univariate_regression_model.hpp"
#include <gsl/gsl_sf_gamma.h>
#include <math.h>

ur_model::ur_model(double alpha,double gamma, double vconst, Data<double> *data)
:probability_model(),m_alpha(alpha),m_gamma(gamma),m_v(vconst)
{
    m_data_ur = data;
    m_data_points = m_data_ur->get_cols();
    m_likelihood_term = -0.5*log(m_v)+m_alpha*log(m_gamma)-gsl_sf_lngamma(m_alpha);
    m_inv_v = 1.0 / m_v;
    m_ysum = new double[m_data_points];
    m_ysum2 = new double[m_data_points];

    m_ysum[0]= m_data_ur->m_X[0][0];
    m_ysum2[0] = (m_data_ur->m_X[0][0])*(m_data_ur->m_X[0][0]);

    for (unsigned int i=1; i<m_data_points; i++){
        m_ysum[i]=m_data_ur->m_X[0][i] + m_ysum[i-1];
        m_ysum2[i]=(m_data_ur->m_X[0][i])*(m_data_ur->m_X[0][i]) +m_ysum2[i-1];
    }
        
  
}

ur_model::~ur_model(){
  //    if(m_data_ur)
  //        delete m_data_ur;
    delete [] m_ysum;
    delete[] m_ysum2;
}
    

double ur_model::log_likelihood_interval(changepoint *obj1, changepoint *obj2, changepoint *objl1){

    double like, y, y2,r;
    unsigned long long int dataindex1, dataindex2;

    like=m_likelihood_term;

  
    //number of observations in each interval
    dataindex1=obj1->getdataindex();
    dataindex2=obj2->getdataindex();
    r = dataindex2-dataindex1;

     if (r<0){
        cout<<"r can not be less than 0"<<endl;
        exit(1);
    }

     if(r==0){    

         like = 0;
     }
     else{

         if(dataindex1==0){
             y=m_ysum[dataindex2-1];
             y2=m_ysum2[dataindex2-1];
         }
         else{
             
             y=m_ysum[dataindex2-1]-m_ysum[dataindex1-1];
             y2=m_ysum2[dataindex2-1]-m_ysum2[dataindex1-1];
         }
	 double r_foo = r / 2.0;
         like += -0.5*log(m_inv_v+r)-(r_foo)*LOG_PI + gsl_sf_lngamma(r_foo+m_alpha)-(r_foo+m_alpha)*log(m_gamma+0.5*(y2-y*y*(1.0/(m_inv_v+r))));
     }

    return(like);

}


double ur_model::calculate_mean(changepoint *obj1, changepoint *obj2, changepoint *objl1){
    
    unsigned long long int dataindex1, dataindex2;

    double mean, r;
    
//number of uncensored observations in each interval

    dataindex1=obj1->getdataindex();
    dataindex2=obj2->getdataindex();
    r = dataindex2-dataindex1;
    
    if (r<0){
        cout<<"r can not be less than 0"<<endl;
        exit(1);
    }

    if (r==0){
        mean=0;
    }
    else{
        
        if(dataindex1==0){
             mean = (1.0/(m_inv_v+r))*(m_ysum[dataindex2-1]);
        }
        else{
            mean = (1.0/(m_inv_v+r))*(m_ysum[dataindex2-1]-m_ysum[dataindex1-1]);
        }
    }
    //    if(r>0)
    //      mean = (m_ysum[dataindex2-1]-m_ysum[dataindex1-1])/double(r);
    return(mean);

}




    
