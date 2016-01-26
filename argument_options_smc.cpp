#include "argument_options_smc.hpp"
#include <time.h>

#include <errno.h>


struct option ArgumentOptionsSMC::lopts[] = {
    {"help", no_argument, NULL, 'h'},
    {"intervals", required_argument, NULL, 'i'},
    {"particles", required_argument, NULL, 'p'},
    {"burnin", required_argument, NULL, 'd'},
    {"thinning", required_argument, NULL, 't'},
    {"movewidth", required_argument, NULL, 'm'},
    {"cpprior", required_argument, NULL, 'n'},
    {"modelprior1", required_argument, NULL, 'a'},
    {"modelprior2", required_argument, NULL, 'b'},
    {"modelpriorur", required_argument, NULL, 'A'},
    {"seed", required_argument, NULL, 's'},
    {"mean", no_argument, NULL, 'l'},
    {"grid", required_argument, NULL, 'g'},
    {"emptyintervals", no_argument, NULL, 'e'},
    {"writecps", no_argument, NULL, 'v'},
    {"model", required_argument, NULL, 'c'},
    {"essthreshold",required_argument,NULL,'f'},
    {"writeess",no_argument,NULL,'w'},
    {"importsampling", no_argument, NULL, 'z'},
    {"spacingprior", no_argument, NULL,  'P'},
    {"rejectionsampling", no_argument, NULL, 'r'},
    {"sequentialmcmc", no_argument, NULL, 'S'},
    {"priorproposal", no_argument, NULL, 'o'},
    {NULL, 0, NULL, 0}
};


ArgumentOptionsSMC::ArgumentOptionsSMC(){
  m_particles = 10000;
  m_sample_sizes = NULL;
  m_num_intervals = 112;
  m_burnin = 1000;
  m_thinning = 10; 
  m_move_width = 0;
  m_cp_prior = 2/(double)112;
  m_gamma_prior_1 = 0.1;
  m_gamma_prior_2 = 0.1;
  m_v = 10;
  srand (time(NULL));
  m_seed = rand() % 10000;
  m_calculate_filtering_mean = 1;
  m_grid = 0;
  m_disallow_empty_intervals_between_cps = 0;
  m_write_cps_to_file = 0;
  m_model = "poisson";
  m_ESS_threshold = 0.5;
  m_print_ESS = 0;
  m_importance_sampling = 0;
  m_rejection_sampling = 0;
  m_spacing_prior = 0;
  m_smcmc=false;
  m_prior_proposals=false;
  m_start=0;
  m_end=112;
  m_datafile="coal_data_renormalised.txt";
}

void ArgumentOptionsSMC::parse(int argc, char * argv[]){

   const char *sopts="hi:p:d:t:m:n:a:b:A:s:lg:evwc:f:zPrSo";

  //Parse arguments
  char opt;
  while ((opt = getopt_long(argc, argv, sopts, lopts, NULL)) != -1) {
    switch(opt){
    case 'i': 
      m_num_intervals = stringtolong(optarg,opt);
      break;
    case 'p':
      m_particles = stringtolong(optarg,opt);
      break;
    case 'd':
      m_burnin = stringtolong(optarg,opt);
      break;
    case 't':
      m_thinning = stringtolong(optarg,opt);
      break;
    case 'm':
      m_move_width = stringtodouble(optarg,opt);
      break;
    case 'n':
      m_cp_prior = stringtodouble(optarg,opt);
      break;
    case 'a':
      m_gamma_prior_1 = stringtodouble(optarg,opt);
      break;
    case 'b':
      m_gamma_prior_2 = stringtodouble(optarg,opt);
      break;
    case 'A':
      m_v = stringtodouble(optarg,opt);
      break;
    case 's':
      m_seed = stringtolong(optarg,opt);
      break;
    case 'l':
      m_calculate_filtering_mean = 1;
      break;
    case 'g':
      m_grid = stringtolong(optarg,opt);
      break;
    case 'e':
      m_disallow_empty_intervals_between_cps = 1;
      break;
    case 'v':
      m_write_cps_to_file = 1;
      break;
    case 'w':
      m_print_ESS = 1;
      break;
    case 'c':
      m_model = optarg;
      break;
    case 'f':
      m_ESS_threshold = stringtodouble(optarg,opt);
      break;
    case 'z':
      m_importance_sampling = 1;
      break;
    case 'r':
      m_rejection_sampling = 1;
      break;
    case 'P':
      m_spacing_prior = 1;
      break;
    case 'S':
      m_smcmc = true;
      break;
    case 'o':
      m_prior_proposals = true;
      break;
    case 'h':
      usage(0,argv[0]);
      break;
    default:
      usage(1,argv[0]);
   
    }
  }

  if(argc-optind != 3){
    usage(1,argv[0]);
  }

  m_datafile = argv[optind++]; 
  m_start = atof(argv[optind++]); 
  m_end = atof(argv[optind++]);

  if(m_grid == 0){
    m_grid = m_num_intervals;
  }else if(!(m_grid >= m_num_intervals and (m_grid % m_num_intervals) == 0)){
    cerr << "the grid has to be multiple of the number of intervals setting to equal the number of intervals" << endl;
    m_grid = m_num_intervals;
  }
  
  if(m_move_width == 0){
    double change_in_time = (m_end-m_start)/(double)m_num_intervals;
    m_move_width = change_in_time/3.0;
  }
}


void ArgumentOptionsSMC::usage(int status,char * programname){
  cerr << endl;
  cerr << "Usage: " << programname <<  " [OPTIONS] DATAFILE STARTTIME ENDTIME" << endl;
  cerr << endl;
  cerr << "-h | --help              print this text and exit" << endl;
  cerr << "-i | --intervals         number of sequential time intervals (default = " << m_num_intervals << ")" << endl;
  cerr << "-p | --particles         number of particles (default = " << m_particles << ")" << endl;
  cerr << "-f | --essthreshold      the threshold at which to resample the particles as a percentage (default = " << m_ESS_threshold << ")" << endl;
  cerr << "-c | --model             can either be sncp (shot noise cox process), poisson (poisson process)," << endl;
  cerr << "                         pregression (poisson regression) or ur (univariate regression) (default = " << m_model << ")" << endl;
  cerr << "-n | --cpprior           the Poisson process prior parameter for the changepoints (default = " << m_cp_prior << ")" << endl;
  cerr << "-a | --modelprior1       prior parameter 1 (dependent on model), see documentation (default = " << m_gamma_prior_1 << ")" << endl;
  cerr << "-b | --modelprior2       prior parameter 2 (dependent on model), see documentation (default = " << m_gamma_prior_2 << ")" << endl;
  cerr << "-A | --modelpriorur      regressor prior variance parameter for univariate regression (default = " << m_v << ")" << endl;
  cerr << "-s | --seed              set the seed for generating random variables (default = current time)" << endl;
  cerr << "-l | --mean              calculate the filtering estimate of the mean over --grid," << endl;
  cerr << "                         no argument required (default = " << m_calculate_filtering_mean << ")" << endl;
  cerr << "-g | --grid              the number of grid points over which to calculate the mean" << endl;
  cerr << "                         must be a multiple of intervals (default = END)" << endl;
  cerr << "-v | --writecps          write changepoints, intensity, and weights at the final time point," << endl;
  cerr << "                         no argument required (default = " << m_write_cps_to_file << ")" << endl;
  cerr << "-w | --writeess          write ESS to file (default = " << m_print_ESS << ")" << endl;
  cerr << "-z | --importsampling    do importance sampling for the coal data, no argument required (default = " << m_importance_sampling << ")" << endl;
  cerr << "-P | --spacingprior      use spacing prior for the coal data, no argument required (default = " << m_spacing_prior << ")" << endl;
  cerr << "-r | --rejectionsampling use rejection sampling rather than mcmc for the coal data," << endl;
  cerr << "                         no argument required (default = " << m_rejection_sampling << ")" << endl;
  cerr << "-S | --sequentialmcmc    use full MCMC on the whole interval [0,t_i] at each update (default = " << m_smcmc << ")" << endl;
  cerr << "-o | --priorproposal     make proposals from prior in each SMC update interval (default = " << m_prior_proposals << ")" << endl;
 
  cerr << endl;
  cerr << "Optional parameters to set when using RJ to sample on each interval" << endl;
  cerr << endl;

  cerr << "-b | --burnin          number of burnin iterations (default = " << m_burnin << ")" << endl;
  cerr << "-t | --thinning        number of iterations to discard between samples (default = " << m_thinning << ")" <<  endl;
  cerr << "-m | --movewidth       the allowable move width on either side of the changepoint (default = (intervalwidth)/20)" << endl;
  cerr << "-e | --emptyintervals  do not allow intervals between changepoints with no datapoints, no argument required (default = " <<  m_disallow_empty_intervals_between_cps << ")" << endl;

 
  cerr << endl;

  cerr << "Example: Poisson process, coal data" << endl;
  cerr << programname << " --intervals 112 --grid 112 --cpprior $(echo '2.0/112.0' | bc -l) --essthreshold 0.3 --writeess --mean -z coal_data_renormalised.txt 0 112" << endl;
  cerr << endl;

  cerr << "Example: shot noise cox process model" << endl;
  cerr << programname << " --model sncp --intervals 40 --particles 500 --movewidth 5 --cpprior .0025 --modelprior1 $(echo '2.0/3.0' | bc -l) --modelprior2 .01 --grid 1000 --essthreshold 0.4  --writeess -z  shot_noise.txt 0 2000" << endl;
  cerr << endl;

  cerr << "Example: Poisson process, simulated data (first run the script simulate_poisson_process.R)" << endl;
  cerr << programname << " --intervals 100 --grid 100 --cpprior 0.0001 --essthreshold 0.3 --writeess --mean -a 1 -b 0.01 -s 0 pp.txt 0 10000" << endl;
  cerr << endl;

  cerr << "Example: Univariate regression (first run the script simulate_gaussian_regression.R)" << endl;
  cerr << programname << " --model ur --intervals 100 --grid 100 --cpprior $(echo '1.0/1000.0' | bc -l) --essthreshold 0.3 --writeess --mean -a 2 -b 2 ur.txt 0 1000" << endl;
  cerr << endl;

  cerr << "Example: Poisson regression  (first run the script simulate_poisson_regression.R)" << endl;
  cerr << programname << " --model ur --intervals 100 --grid 100 --cpprior $(echo '1.0/1000.0' | bc -l) --essthreshold 0.3 --writeess --mean -a .1 -b .1 ur.txt 0 1000" << endl;
  exit(status);

}

double ArgumentOptionsSMC::stringtodouble(char * parameter, char opt){
  char * pend;
  errno=0;
  double foo = strtod(parameter,&pend); 
  errorchecking(pend,opt); 
  return foo;
}

long int ArgumentOptionsSMC::stringtolong(char * parameter, char opt){
  char * pend;
  errno=0;
  long int foo = strtol(parameter,&pend,10); 
  errorchecking(pend,opt); 
  return foo;
}


void ArgumentOptionsSMC::errorchecking(char * pend, char opt){
  
  if (*pend != '\0' || errno == ERANGE){ 
    int i=0;
    while(lopts[i].name != NULL){
      if(lopts[i].val == opt){
	cerr << endl;
	cerr << "invalid argument for " << lopts[i].name << endl;
	exit(1);
      }
      i++;
    }
  }
}
		       
