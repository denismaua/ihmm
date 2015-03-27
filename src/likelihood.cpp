/* 
 * imprecise hidden Markov model
 *
 * upper and lower log-likelihood inference algorithm
 * 
 * Compute the upper and lower probabilities of a
 * separately specified imprecise hidden Markov model
 * for a given sequence of observations
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
#include <utility>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <unistd.h>

#include "inference.h"
#include "utils.h"
#include "version.h"

void usage(char*); // print usage message

int main(int argc, char* argv[]) {

  /* Print program info */
  print_info("Log-likelihood Inference"); // show distribution info

  cout.precision(15); // set high output precision    
	
  /* Model specification */
  unsigned int N = 0, M = 0;
  NUM *Al, *Au, *Bl, *Bu, *Pl, *Pu;
  /* Input length */
  unsigned int K = 0;
  /* observation sequence */
  unsigned int *O;
	
  /* Parse command line with getopt */
  int c;
  opterr = 0;
  int debug = 0; // output additional information
  int opt = 1; // use quick select as optimization routine
  int usesecs = 0; // report time in seconds
	
  while ((c = getopt (argc, argv, "dho:s")) != EOF)
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
	
  char *ifilename = argv[optind] ;
  char *dfilename = argv[optind+1];
  char *ofilename = argv[optind+2];
	
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
  cout << endl << endl;
	
  // load model from file
  bool error = load_model(Al, Au, Bl, Bu, Pl, Pu, N, M, ifilename);
  if (error) {
    cout << "Error loading model from file" << endl;
    return 1;
  }

  cout << "N: " << N << endl;
  cout << "M: " << M << endl;
	
  // read observations from file
    
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
    // run algorithm
    pair<NUM,NUM> ll = likelihood(Al,Au,Bl,Bu,Pl,Pu,O,N,M,K,debug);
    iend = getcputime();         // record finishing time
        
    // output result
    cout << "LL: [" << ll.first << ", " << ll.second << "]" << endl;
    ofs << ll.first << " " << ll.second << endl;


    // Compute inference running time
    cout << "Inference running time: ";
    if (usesecs)  // report time in secs
      cout << fixed << (iend-istart) << "s" <<  endl;
    else { // format time
      difftime(istart,iend,time);
      cout << fixed << time[0] << "h " << time[1] << "m " << time[2] << "s" <<  endl;
    }
    cout << endl;
        
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
  cout << "Usage: " << program_name << " [-dho:s] model_filename data_filename output_filename\n\n";
  cout << "\t-d\t debug mode               [off]\n";
  cout << "\t-h\t Print out usage information (this)\n";
  //cout << "\t-n\t do not use scaling      [off]\n"; // not implemented
  cout << "\t-o\t optimization algorithm  [1]\n";
  cout << "\t-s\t report time in secs     [off]\n";
  cout << endl;
}
