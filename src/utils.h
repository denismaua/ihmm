/* 
 * imprecise hidden Markov model
 *
 * Utility functions
 * 
 * See COPYRIGHT and VERSION files for copyright and version info
 *
 */

#ifndef _UTILS_H
#define _UTILS_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include "inference.h"

typedef std::vector< unsigned int > sequence; // data sequence type
// input
bool load_model(NUM*&, NUM*&, NUM*&, NUM*&, NUM*&, NUM*&, unsigned int&, unsigned int&, char*); // load uniform model from file
bool load_model2(NUM**&, NUM**&, NUM**&, NUM**&, NUM*&, NUM*&, unsigned int* &, unsigned int* &, unsigned int&, char*); // load non-uniform model from file
bool load_data(std::vector<sequence>&, unsigned int, unsigned int, char*); // load data from file
bool load_data2(std::vector<sequence>&, unsigned int, unsigned int*, char*); // load data from file (non-uniform model)
// output
void print_precise_probabilities (NUM **, NUM **, NUM *, unsigned int*, unsigned int*, unsigned int);  // print probabilities
void print_probabilities (NUM*, NUM*, NUM*, NUM*, NUM*, NUM*, unsigned int, unsigned int); // print credal sets (uniform)
void print_probabilities2 (NUM**, NUM**, NUM**, NUM**, NUM*, NUM*, unsigned int*, unsigned int*, unsigned int); // print credal sets (non-uniform)
void print_observations (std::vector<sequence>); // print observations
// time measurement
double getcputime(void); // get cpu time for performance measurement
void difftime(double, double, unsigned int*); // format time
// validation
bool is_reachable(NUM*, NUM*, NUM); // check reachability
#endif
