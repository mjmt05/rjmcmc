#include "argument_options_vastdata.hpp"
#include <time.h>

#include <errno.h>


struct option ArgumentOptionsVast::lopts[] = {
    {"help", no_argument, NULL, 'h'},
    {"intervals", required_argument, NULL, 'i'},
    {"particles", required_argument, NULL, 'p'},
    {"burnin", required_argument, NULL, 'd'},
    {"thinning", required_argument, NULL, 't'},
    {"movewidth", required_argument, NULL, 'm'},
    {"cpprior", required_argument, NULL, 'n'},
    {"modelprior1", required_argument, NULL, 'a'},
    {"modelprior2", required_argument, NULL, 'b'},
    {"seed", required_argument, NULL, 's'},
    {"mean", no_argument, NULL, 'l'},
    {"grid", required_argument, NULL, 'g'},
    {"emptyintervals", no_argument, NULL, 'e'},
    {"essthreshold",required_argument,NULL,'f'},
    {"writeess",no_argument,NULL,'w'},
    {"num_bins",required_argument,NULL,'B'},
    {"loss_function", required_argument, NULL, 'L'},
    {"min_iterations", required_argument, NULL, 'M'},
    {"fixed_sample", no_argument, NULL, 'F'},
    {"sample_sizes", required_argument, NULL, 'S'},
    {NULL, 0, NULL, 0}
};


ArgumentOptionsVast::ArgumentOptionsVast(){
  m_particles = 10000;
  m_num_intervals = 240;
  m_burnin = 5000;
  m_thinning = 10; 
  m_move_width = 0.01;
  m_cp_prior = 0.005;
  m_gamma_prior_1 = 0.05;
  m_gamma_prior_2 = 0.05;
  srand (time(NULL));
  m_seed = rand() % 10000;
  m_calculate_filtering_mean = 1;
  m_grid = 0;
  m_disallow_empty_intervals_between_cps = 0;
  m_ESS_threshold = 0.5;
  m_print_ESS = 0;
  m_num_bins = 50;
  m_loss_type = AVERAGE;
  m_min_iterations = 500;
  m_fixed_sample_size = false;
  m_sample_sizes = "";
 }

void ArgumentOptionsVast::parse(int argc, char * argv[]){

   const char *sopts="hi:p:d:t:m:n:a:b:s:lg:evwf:B:L:M:F";

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
    case 'w':
      m_print_ESS = 1;
      break;
    case 'f':
      m_ESS_threshold = stringtodouble(optarg,opt);
      break;
    case 'F':
      m_fixed_sample_size = true;
      break;
    case 'h':
      usage(0,argv[0]);
      break;
    case 'B':
      m_num_bins = stringtolong(optarg,opt);
      break;
    case 'L':
      m_loss_type = static_cast<Loss_Function>(atoi(optarg));
      break;
    case 'M':
      m_min_iterations = stringtolong(optarg, opt);
      break;
    case 'S':
      m_sample_sizes = optarg;
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


void ArgumentOptionsVast::usage(int status,char * programname){
  cerr << endl;
  cerr << "Usage: " << programname <<  " [OPTIONS] DATAFILE STARTTIME ENDTIME" << endl;
  cerr << "Note: DATAFILE should be a file with the list of files (including path) for all the processes" << endl;
  cerr << endl;
  cerr << "-h | --help              print this text and exit" << endl;
  cerr << "-F | --fixed_sample      do fixed sample size rather than adaptive sample size, no argument required (default = " << m_fixed_sample_size << ")" << endl;
  cerr << "-i | --intervals         number of sequential time intervals (default = " << m_num_intervals << ")" << endl;
  cerr << "-p | --particles         number of particles. For adaptive sample size the total number of samples will be " << endl;
  cerr << "                         particles multiplied by the number of processes (default = " << m_particles << ")" << endl;
  cerr << "-f | --essthreshold      the threshold at which to resample the particles as a percentage (default = " << m_ESS_threshold << ")" << endl;
  cerr << "-n | --cpprior           the Poisson process prior parameter for the changepoints (default = " << m_cp_prior << ")" << endl;
  cerr << "-a | --modelprior1       prior parameter 1 for Poisson process, see documentation (default = " << m_gamma_prior_1 << ")" << endl;
  cerr << "-b | --modelprior2       prior parameter 2 for Poisson process, see documentation (default = " << m_gamma_prior_2 << ")" << endl;
  cerr << "-s | --seed              set the seed for generating random variables (default = current time)" << endl;
  cerr << "-l | --mean              calculate the filtering estimate of the mean over --grid," << endl;
  cerr << "                         no argument required (default = " << m_calculate_filtering_mean << ")" << endl;
  cerr << "-g | --grid              the number of grid points over which to calculate the mean" << endl;
  cerr << "                         must be a multiple of intervals (default = END)" << endl;
  cerr << "-w | --writeess          write ESS to file, no argument required (default = " << m_print_ESS << ")" << endl;

  cerr << endl;

  cerr << "Optional parameters for adaptive sample size" << endl;
  cerr << endl;

  cerr << "-M | --min_iterations  minimum number of samples to allocate to each individual on each interval (default = " << m_min_iterations << ")" << endl;
  cerr << "-L | --loss_function   0 for minimax loss or 1 for average loss (default = " << m_loss_type << ")" << endl;
  cerr << "-B | --num_bins        number of bins to use when constructing a histogram of sampled changepoints (default = " << m_num_bins << ")" << endl;
  cerr << "-S | --sample_sizes    a file with the a set of sample sizes to be read in for an individual" << endl;

  cerr << endl;
  cerr << "Optional parameters to set when using RJ to sample on each interval" << endl;
  cerr << endl;

  cerr << "-b | --burnin          number of burnin iterations (default = " << m_burnin << ")" << endl;
  cerr << "-t | --thinning        number of iterations to discard between samples (default = " << m_thinning << ")" <<  endl;
  cerr << "-m | --movewidth       the allowable move width on either side of the changepoint (default = " << m_move_width << ")" << endl;
  cerr << "-e | --emptyintervals  do not allow intervals between changepoints with no datapoints, no argument required (default = " <<  m_disallow_empty_intervals_between_cps << ")" << endl;

 
  cerr << endl;

  cerr << "Example: Vast Data" << endl;
  cerr << programname << " filenames.txt 0 10" << endl;
  cerr << endl;

  exit(status);

}

double ArgumentOptionsVast::stringtodouble(char * parameter, char opt){
  char * pend;
  errno=0;
  double foo = strtod(parameter,&pend); 
  errorchecking(pend,opt); 
  return foo;
}

long int ArgumentOptionsVast::stringtolong(char * parameter, char opt){
  char * pend;
  errno=0;
  long int foo = strtol(parameter,&pend,10); 
  errorchecking(pend,opt); 
  return foo;
}


void ArgumentOptionsVast::errorchecking(char * pend, char opt){
  
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
		       
