#ifndef SMC_PP_MCMC_NC_HPP
#define SMC_PP_MCMC_NC_HPP


#include "changepoint.hpp"
#include "SMC_PP.hpp"
#include "RJMCMC_PP.hpp"
#include "probability_model.hpp"
#include "function_of_interest.hpp"
#include "rejection_sampling.hpp"
#include <list>
#include <utility>  

class SMC_PP_MCMC : public SMC_PP<changepoint>{
public:

  SMC_PP_MCMC(double = 0, double = 1, unsigned int = 1, int =0, int=0, double=1, double=0, probability_model ** =0,int=1, bool=0,bool=0, bool=0, bool=0, bool = 0, int=0);
    ~SMC_PP_MCMC();
    virtual void sample_particles(double, double);
    virtual void resample_particles(double,double,int, const char *,int);
    virtual void ESS_resample_particles(double,int);
    virtual void calculate_function_of_interest(double, double);
    virtual void calculate_weights_join_particles(int,int);
    virtual void delete_samples(int);

    void print_exp(const char *);
    void print_var_exp(int,const char *);
    void print_intensity(int=0,const char * = "intensity.txt");
    void print_prob(const char *);
    void print_exp_sequential(int, const char *);
    void print_prob_sequential(int, const char *);
    void set_look_back(bool lb);
    void increase_A_particles(int,unsigned long long int, unsigned long long int *);
  
    void store_sample_sizes();
    void print_variable_sample_sizes(const char *);
    /*always initialise function of interest before setting look back*/
    void initialise_function_of_interest(int, bool, bool, bool=0, double=0, bool=0);
    unsigned long long int ** get_sample_sizges(){return m_sample_sizes;};
    unsigned long long int * get_min_sample_sizes(){return m_min_sample_size;};
    Function_of_Interest * get_function_of_interest(int ds){return m_functionofinterest[ds];}
    long double * get_prob(int ds){return m_functionofinterest[ds]->get_prob();}
    long double * get_variance_g(int ds){return m_functionofinterest[ds]->get_variance_g();}
    long double * get_g(int ds){return m_functionofinterest[ds]->get_g();}
    void set_neighbouring_intervals(bool ni){m_neighbouring_intervals=ni;}
    void set_RJ_parameters(int, int, double, const char *proposaltype=NULL, void * v=NULL);
    void set_variable_parameters(Divergence_Type dt, Loss_Function lt, int nb,int ii, long long int mll, int fg=0){m_divergence_type=dt; m_loss_type=lt, m_num_of_bins=nb; m_initial_iterations=ii; m_max_lookup_length=mll;m_foi_grid=fg;}
    void set_proposal_prior(double pp){m_proposal_prior=pp; m_prior_diff=log(m_nu)-log(pp);}
    void non_conjugate(){m_conjugate=0;}
    void use_spacing_prior(double space=0);
    void print_zero_weights(int, const char *);
    void print_rejection_sampling_acceptance_rates(int, const char *);
  void sample_intensities(Particle<changepoint> **, double, unsigned int, int);
  void set_discrete_model(){m_discrete = true;}
  

private:
       
        rj_pp ** m_rj_B;
        rejection_sampling ** m_rejection_sampling;
        rj_pp ** m_rj_A;
	int m_length_grid;
	bool m_calculate_intensity;
        bool m_do_exact_sampling;
        bool m_use_spacing_prior;
	Function_of_Interest ** m_functionofinterest;
	bool m_do_SMC_past;
	int m_initial_iterations;
	double * m_vec_KLS;
	bool m_store_sample_sizes;
	unsigned long long int ** m_sample_sizes;
	unsigned long long int * m_min_sample_size;
        unsigned long long int * m_current_sample_size;
	const char *m_proposal_type;
	void * m_vec_proposal_type;
	bool m_neighbouring_intervals;
	double m_nu, m_var_nu;
       //parameter for enforced spacing when doing rejection sampling;
        double m_spacing_prior;
        double m_cp_start;
        
	//RJ parameters when sampling on the intervals
	int m_thin;
	double m_move_width;
	int m_burnin;
	double	m_prior_diff;
	double m_proposal_prior;
	//variable m parameters
	int m_num_of_bins;
	Divergence_Type m_divergence_type;
	Loss_Function m_loss_type;
	int m_foi_grid;
	unsigned long long int m_max_lookup_length;
	bool m_conjugate;
        double **m_rejection_sampling_acceptance_rate;
        unsigned int **m_num_zero_weights;
  bool m_discrete;
       
 
  void increase_vector(int, unsigned long long int);
	static bool MyDataSort(const pair<double,int>&, const pair<double,int>&);

};




#endif
