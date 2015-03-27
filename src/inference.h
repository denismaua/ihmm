/* 
 * imprecise hidden Markov model
 *
 * Inference algorithms
 * 
 * See COPYRIGHT and VERSION files for copyright and version info
 *
 */

#ifndef _INFERENCE_H
#define _INFERENCE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

typedef double NUM; // data type for real numbers

// Data structs
struct associate
{
  /* real, integer pair
     useful for  bookeeping vector position information
     when ordering by value
  */
    NUM value;
    int index;
    associate (NUM v, unsigned int i) : value(v), index(i) {};
    associate () : value((NUM)0), index(0) {};
};

inline bool operator<(const associate &x, const associate &y)
{
  return x.value < y.value;
}
inline bool operator<= (const associate &x, const associate &y)
{
  return x.value <= y.value;
}
inline bool operator>(const associate &x, const associate &y)
{
  return x.value > y.value;
}
inline bool operator>=(const associate &x, const associate &y)
{
  return x.value >= y.value;
}
inline bool operator==(const associate &x, const associate &y)
{
  return x.value == y.value;
}
inline bool operator!=(const associate &x, const associate &y)
{
  return x.value != y.value;
}
inline std::ostream& operator<< (std::ostream &s, const associate &x)
{
    return s << '<' << x.value << ';' << x.index << '>';
}

// Methods

unsigned int median( std::vector<associate, std::allocator<associate> >&, unsigned int, unsigned int, unsigned int); // median of three
NUM knapsack(NUM*, NUM*, NUM*, unsigned int, unsigned int, NUM* x=NULL); // continuous knapsack
NUM greedyknapsack(NUM*, NUM*, NUM*, unsigned int, NUM* x=NULL); // greedy knapsack solver
NUM heapknapsack(NUM*, NUM*, NUM*, unsigned int, NUM* x=NULL); // greedy knapsack with partial heapsort
NUM qsknapsack(NUM*, NUM*, NUM*, unsigned int, NUM* x=NULL); // quick select knapsack solver
NUM quickselectknapsack( std::vector<associate, std::allocator<associate> > &, NUM*, unsigned int, unsigned int, NUM, unsigned int, NUM* x=NULL); 
NUM viterbi(NUM*, NUM*, NUM*, unsigned int*, unsigned int, unsigned int, unsigned int, unsigned int*, int=0, int=1); // most probable explanation (viterbi) [homogeneous model]
NUM viterbi2(NUM**, NUM**, NUM*, unsigned int*, unsigned int*, unsigned int*, unsigned int, unsigned int*, int=0, int=1); // most likely state sequence (viterbi) [non-homogeneous model]
std::pair<NUM,NUM> likelihood(NUM*, NUM*, NUM*, NUM*, NUM*, NUM*, unsigned int*, unsigned int, unsigned int, unsigned int, int=0, int=1); // likelihood interval [homogeneous models]
std::pair<NUM,NUM> likelihood2(NUM**, NUM**, NUM**, NUM**, NUM*, NUM*, unsigned int*, unsigned int*, unsigned int*, unsigned int, int=0, int=1); // likelihood interval [non-homogenous models]
std::pair<NUM,NUM> filtering(unsigned, NUM*, NUM*, NUM*, NUM*, NUM*, NUM*, unsigned int*, unsigned int, unsigned int, unsigned int, double, int=0, int=1); // predictive probability interval [homogeneous models]
std::pair<NUM,NUM> filtering2(unsigned, NUM**, NUM**, NUM**, NUM**, NUM*, NUM*, unsigned int*, unsigned int*, unsigned int*, unsigned int, double, int=0, int=1); // predictive probability interval [non-homogeneous models]

#endif
