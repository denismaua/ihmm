/* 
 * imprecise hidden Markov model
 *
 * predictive inference (filtering)
 * 
 * Computes the lower an upper probabilities of the current hidden state
 * for a given sequence of observations.
 *
 * non-homogeneous models.
 *
 * See COPYRIGHT and VERSION files for copyright and version info
 *
 */
using namespace std;
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <unistd.h>


#include "inference.h"
#include "version.h"
#include "utils.h"

void usage(char*); // print usage message

int main(int argc, char* argv[]) {

  /* Print program info */
  print_info("Filtering"); // show distribution info

  cout.precision(15); // set high output precision
	    
	
  /* Model specification */
  unsigned int* N; unsigned int* M; // number os states and symbols (per step)
  NUM **Al, **Au, **Bl, **Bu; // interval probability matrices
  NUM *Pl, *Pu; 
  /* Input length */
  unsigned int K=0;
  /* observation sequence */
  unsigned int *O;
    
  double epsilon = 1E-6; // binary search resolution
	

  /* Parse command line with getopt */
  int c;
  opterr = 0;

  // default parameter values
  int debug = 0; // print out extra info?
  int opt = 1; // use quick select as optimization routine
  int usesecs = 0; // report time in seconds

  while ((c = getopt (argc, argv, "de:ho:s")) != -1)
    switch (c)
      {
      case 'd':
	// debug mode
	debug = 1;
	break;
      case 'e':
	epsilon = atof(optarg);
	break;
      case 'h':
	usage(argv[0]);
	return 1;
	break;		
      case 'o':
	opt = atoi(optarg);
	if (opt > 2) opt = 1;
	if (opt < 0) opt = 0;
	break;
      case 's':
	// report in in secs
	usesecs = 1;
	break;
      case '?':
	if (optopt == 'e')
	  fprintf (stderr, "Option -%c requires a double value as argument.\n", optopt);
	else if (optopt == 'o')
	  fprintf (stderr, "Option -%c requires an integer value as argument.\n", optopt);
	else if (isprint (optopt))
	  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
	  fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	return 1;
      default:
	abort ();
      }	
	
  if (argc-optind < 3) {
    usage(argv[0]);
    return 1;
  }
	
  char *ifilename = argv[optind];
  char *dfilename = argv[optind+1];
  char *ofilename = argv[optind+2];

  // report options
  cout << "model filename: " << ifilename << endl;
  cout << " data filename: " << dfilename << endl;
  cout << "        solver: ";
  if (opt == 0)
    cout << "greedy knapsack";        
  else if (opt == 1)
    cout << "quick select knapsack";
  else if (opt == 2)
    cout << "greedy heapsort knapsack";
  else
    cout << "unknown";
  cout << endl;
  cout << endl << endl;
	
  // load model from file
  bool error = load_model2(Al, Au, Bl, Bu, Pl, Pu, N, M, K, ifilename);
  if (error) {
    cout << "Error loading model from file" << endl;
    return 1;
  }
  
  cout << "T: " << K << endl;
  cout << "N: "; for (unsigned n=0;n<K;n++) cout << N[n] << " ";
  cout << endl;
  cout << "M: "; for (unsigned n=0;n<K;n++) cout << M[n] << " ";
  cout  << endl;
	
  vector<sequence> D; // dataset
  // read observations from file
  error = load_data2(D, K, M, dfilename);
  if (error) {
    cout << "Error loading data from file" << endl;
    return 1;
  }

  cout << "D: " << D.size() << endl << endl; // nr. of sequences read
	
  if (debug) {
		
    print_probabilities2(Al, Au, Bl, Bu, Pl, Pu, N, M, K);
    print_observations(D);

  } 
  cout << endl;
  
  /* for elapsed time measuring */
  unsigned int time[3];	
  double istart, iend;
	
  ofstream ofs; // open stream
  ofs.open(ofilename);
  if ( !ofs ) { // file couldn't be opened
    cerr << "Error: file " << ofilename << " could not be opened" << endl;
    return 1;
  }
  
  for (unsigned int s=0; s<D.size(); ++s) {
       
    K = D[s].size();
    O = new unsigned int[K];
    if (O==NULL) {
      cerr << "Could not initialize observation vector!" << endl;
      return 1;	
    }
        
    unsigned int t=0;
    for (sequence::iterator o = D[s].begin(); o != D[s].end(); ++o,++t)
      O[t] = *o;

    cout << "#" << (s+1) <<"/" << D.size() << endl;
    cout << "T: " << K << endl;
	
    istart = getcputime();       // record starting time
    for (unsigned q=0; q<N[K-1]; q++)
      {
     if (debug) std::cout << std::endl << "Q=" << q << std::endl;
    	// run algorithm
    	pair<NUM,NUM> prob = filtering2(q,Al,Au,Bl,Bu,Pl,Pu,O,N,M,K,epsilon,debug,opt);
    	// output result
    	cout << scientific << "PROB(Q" << K << "=" << q << "|DATA): [" << prob.first << ", " << prob.second << "]" << endl;
    	ofs << prob.first << " " << prob.second << " ";
      }
    iend = getcputime();

    // Compute inference running time
    cout << "Inference running time: ";
    if (usesecs)  // report time in secs
      cout << fixed << (iend-istart) << "s" <<  endl;
    else { // format time
      difftime(istart,iend,time);
      cout << fixed << time[0] << "h " << time[1] << "m " << time[2] << "s" <<  endl;
    }
    cout << endl;
    ofs << endl;
        
    // free memory
    delete [] O;
    
  }
  ofs.close();
        
	
  // Total elapsed time
  double etime = getcputime();
  difftime(0,etime,time);
    
  cout << "Total running time: " << fixed << time[0] << "h " << time[1] << "m " << time[2] << "s" <<  endl;
  cout << endl;
    
  // free memory

  return 0;
}

void usage(char* program_name) {
  /* Print out usage message */
  cout << "Usage: " << program_name << " [-de:ho:s] model_filename data_filename output_filename\n\n";
  cout << "\t-d\t debug mode              [off]\n";
  cout << "\t-e\t search precision        [0.00000001]\n";
  /* cout << "\t-n\t do not use scaling      [on]\n"; // not implemented */
  cout << "\t-o\t optimization algorithm  [1]\n";
  /* 0: greedy knapsack   1: quickselect knapsack      2: heapsort knapsack */
  cout << "\t-s\t report time in secs     [off]\n";
  cout << "\t-h\t Print out usage information (this)\n";
  cout << endl;
}
