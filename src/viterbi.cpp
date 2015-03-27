/* 
 * imprecise hidden Markov model
 *
 * joint maximin and maximax MPE inference algorithm
 * 
 * Compute the most likely state sequence of a given 
 * separately specified imprecise hidden Markov model
 * for a given sequence of observations using joint
 * probability criteria
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
  print_info("Most Likely State Sequence"); // show distribution info

  cout.precision(15); // set high output precision
	
  /* Model specification */
  unsigned int N, M;
  NUM *Al, *Au, *Bl, *Bu, *Pl, *Pu;
  /* Input length */
  unsigned int K;
  /* observation sequence */
  unsigned int *O;
    
  /* Parse command line with getopt */
  int c;
  opterr = 0;

  // default parameter values
  int debug = 0; // print out extra info?
  int maxmax = 0; // use maxmax mpe instead of maxmin
  int opt = 1; // use quick select as optimization routine
  int usesecs = 0; // report time in seconds
  bool scaling = 1; /* use scaling - not implemented
			0 - no
			1 - yes      */

  while ((c = getopt (argc, argv, "dhmno:s")) != -1)
    switch (c)
      {
      case 'd':
	// debug mode
	debug = 1;
	break;
      case 'h':
	usage(argv[0]);
	return 1;
	break;		
      case 'm':
	// use maxmax instead of maxmin
	maxmax = 1;
	break;
      case 'n':
	// scaling
	scaling = 0;
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
	if (optopt == 'o')
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
  cout << "      criteria: ";
  if (maxmax)
    cout << "maximax";
  else
    cout << "maximin";
  cout << endl;
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
  cout << "       scaling: ";
  if (scaling)
    cout << "on";        
  else if (maxmax)
    cout << "off";
  cout << endl << endl;
	
  // load model from file
  bool error = load_model(Al, Au, Bl, Bu, Pl, Pu, N, M, ifilename);
  if (error) {
    cout << "Error loading model from file" << endl;
    return 1;
  }

  cout << "N: " << N << endl;
  cout << "M: " << M << endl;
	
  vector<sequence> D; // dataset
  // read observations from file
  error = load_data(D, K, M, dfilename);
  if (error) {
    cout << "Error loading data from file" << endl;
    return 1;
  }

  cout << "D: " << D.size() << endl << endl; // nr. of sequences read
	
  if (debug) {
		
    print_probabilities(Al, Au, Bl, Bu, Pl, Pu, N, M);
    print_observations(D);

  } 
  cout << endl;

  if (!is_reachable(Pl,Pu,N))
    cerr << "*** Warning: prior probability intervals are not reachable *** " << endl << endl;

  for (unsigned int i=0; i < N; ++i) {
    NUM* l = Al + i*N; // select local transition credal set A_i
    NUM* u = Au + i*N;
    if (!is_reachable(l,u,N))
      cerr << "*** Warning: transition probability intervals are not reachable ***" << endl << endl;
  }
  
  for (unsigned int i=0; i < N; ++i) {
    // select local emission credal set B_i
    NUM* l = new NUM[M];
    for (unsigned j=0; j<M; j++)
      l[j] = Bl[j*N+i]; 
    NUM* u = new NUM[M];
    for (unsigned j=0; j<M; j++)
      u[j] = Bu[j*N+i]; 
    if (!is_reachable(l,u,M)) {
      cerr << "*** Warning: emission probability intervals are not reachable ***" << endl << endl;
    }
  }
  
  /* for elapsed time measuring */
  unsigned int time[3];	
  double istart, iend;
  NUM prob;
	
  ofstream ofs; // open stream
  ofs.open(ofilename);
  if ( !ofs ) { // file couldn't be opened
    cerr << "Error: file " << ofilename << " could not be opened" << endl;
    return 1;
  }
        
  unsigned int *q; // most likely state sequence

  for (unsigned int s=0; s<D.size(); ++s) {
       
    K = D[s].size();
    q = new unsigned int[K];
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
    if (maxmax) // run viterbi algorithm with upper bounds
      prob = viterbi(Au,Bu,Pu,O,N,M,K,q,debug,opt); 
    else // run viterbi algorithm with lower bounds
      prob = viterbi(Al,Bl,Pl,O,N,M,K,q,debug,opt); 
    iend = getcputime();
    if (!scaling)
      cout << scientific << "PROB: " << exp(prob) << endl;         
    else
      cout << scientific << "logPROB: " << prob << endl;
    // Compute inference running time
    cout << "Inference running time: ";
    if (usesecs)  // report time in secs
      cout << fixed << (iend-istart) << "s" <<  endl;
    else { // format time
      difftime(istart,iend,time);
      cout << fixed << time[0] << "h " << time[1] << "m " << time[2] << "s" <<  endl;
    }
    cout << endl;
    // output result

    // save estimated state sequence
    for (t=0; t < K; ++t)
      ofs << q[t] << " ";
    ofs << endl;
        
    // free memory
    delete [] q;
    delete [] O;
    
  }
  ofs.close();
        
	
  // Total elapsed time
  double etime = getcputime();
  difftime(0,etime,time);
    
  cout << "Total running time: " << fixed << time[0] << "h " << time[1] << "m " << time[2] << "s" <<  endl;
  cout << endl;
    
  // free memory
  delete [] Al;
  delete [] Au;
  delete [] Bl;
  delete [] Bu;
  delete [] Pl;
  delete [] Pu;

  return 0;
}

void usage(char* program_name) {
  /* Print out usage message */
  cout << "Usage: " << program_name << " [-dhmo:s] model_filename data_filename output_filename\n\n";
  cout << "\t-d\t debug mode              [off]\n";
  cout << "\t-m\t use maximax criteria    [off]\n";
  cout << "\t-n\t do not use scaling      [off]\n";
  cout << "\t-o\t optimization algorithm  [1]\n";
  /* 0: greedy knapsack   1: quickselect knapsack      2: heapsort knapsack */
  cout << "\t-s\t report time in secs     [off]\n";
  cout << "\t-h\t Print out usage information (this)\n";
  cout << endl;
}
