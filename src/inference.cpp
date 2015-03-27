/* 
 * imprecise hidden Markov model
 *
 * Inference algorithms
 * 
 * See COPYRIGHT and VERSION files for copyright and version info
 *
 */

#include "inference.h"

unsigned int median( std::vector<associate, std::allocator<associate> > &v, unsigned a, unsigned b, unsigned c) {
  // computes the median of three points at positions v[a], v[b] and v[c]
  if(v[a] <= v[b] && v[b] <= v[c])
    return b;
  if(v[c] <= v[b] && v[b] <= v[a])
    return b;
  if(v[b] <= v[a] && v[a] <= v[c]) 
    return a;
  if(v[c] <= v[a] && v[a] <= v[b]) 
    return a;
  return c;
}

NUM knapsack(NUM l[], NUM u[], NUM w[], unsigned int N, unsigned int method, NUM* x) {
  /* Interface to Continuous Knapsack Solvers
	 
     l; lower bound vector
     u; upper bound vector
     w; value vector
     N; size
     method; solver to use [1-3]
     x; solution (if pointer different from NULL)

  */
  if (method == 1) // use quick select knapsack
    return qsknapsack(l,u,w,N,x);
  if (method == 2) // use heapsort knapsack
    return heapknapsack(l,u,w,N,x);
  else // use greedy knapsack
    return greedyknapsack(l,u,w,N,x);

}

NUM qsknapsack(NUM l[], NUM u[], NUM w[], unsigned int N, NUM* x) {
  /* Remove lower bounds so that QuickSelectKnapsack can be used
	 
     l; lower bound vector
     u; upper bound vector
     w; value vector
     N; size
     x; solution (if pointer different from NULL)

  */
  // map weights to value-index datastruct
  std::vector<associate, std::allocator<associate> > iw;
  NUM C = 0.0; // capacity
  NUM O = 0.0; // optimum
  iw.resize(N);
  NUM* uu = new NUM[N];
  for (unsigned int i = 0; i < N; ++i) {
    iw[i] = associate(w[i],i);
    C += l[i]; // compute lower bound sum
    O += l[i]*w[i]; // fill pack with minimum amounts
    uu[i] = u[i] - l[i]; // fill packs with minimum
    if (x != NULL) x[i] = l[i];
  }

  if (1.0-C > 1e-10) // for precision incosistency
    O = O + quickselectknapsack(iw,uu,0,N,1.0-C,N,x); /* call QuickSelectKnapsack
						                                              with knapasack filled with
											      lower bounds */
  delete [] uu;
  return O;
}

NUM quickselectknapsack( std::vector<associate, std::allocator<associate> > &v, NUM *b, unsigned int n, unsigned int m, NUM C, unsigned int N, NUM* x) {
  /* Random breakpoint search implementation of the continuous knapsack problem (quick select alike)
	 
     runs in expected polynomial time and worst case O(n^2)

     b; upper bound vector
     v; value vector
     n; left element index
     m; right element index + 1
     C; capacity
     N; size
     x; solution (if pointer different from NULL)
  */

  assert(n < m);

  unsigned int i=n+1, j=m-1;
  unsigned int middle = (n+m)/2;
  /* heuristic to choose good pivot (breakpoint)
     Gurwitz, 1992. IEEE TRANSACTIONS ON EDUCATION,
     VOL. 35, NO. 3.
  */
  unsigned int pivot = median(v, n, j, middle);
  associate pValue;

  // move pivot to beggining
  pValue = v[pivot];
  v[pivot] = v[n];
  v[n] = pValue;

  /* rearrange v such that
     v[k] <= v[pivot] for k < pivot
     and v[k] > v[pivot] for k > pivot
  */
  while(i<=j) {
    if(v[i] <= pValue) i++; // v[i] is already in left  half
    else if(v[j] > pValue) j--; // v[j] is already in right half
    else {
      // swap v[i] and v[j]
      associate t;
      t=v[i];
      v[i]=v[j];
      v[j]=t;
      --j;
      ++i;
    }
  }
  v[n]=v[j]; // move pivot to right position
  v[j]=pValue;  

  // j is the breakpoint
  NUM c = 0.0; // sum of weights in the right half
  for(i=j+1; i<m; i++)
    c += b[v[i].index];

  // if c <= C < c+v[j] then we have found a critical point j
  if(c <= C && c+b[v[j].index] > C) {
    /* fill knapsack with maximum for all
       elements in the right half but j
       and assign remaining capacity to j
    */
    NUM res = 0.0;
    for(i=j+1; i<N; i++) {
      if (x != NULL) x[v[i].index] += b[v[i].index]; 
      res += b[v[i].index] * v[i].value;
    }
    if (x != NULL) x[v[j].index] += C-c;
    //std::cout << "j: " << v[j].index << std::endl;
    return res + (C-c)*v[j].value;

  } else {
    // otherwise, keep searching on the correct half
    if(c > C) // breakpoint is too large, search on right half
      return quickselectknapsack(v,b,j+1,m,C,N,x);
    else // breakpoint is too small, search on left half
      return quickselectknapsack(v,b,n,j,C-c-b[v[j].index],N,x);
  }
}

NUM greedyknapsack(NUM l[], NUM u[], NUM w[], unsigned int N, NUM* x) {
  /* Greedy implementation of the continuous knapsack problem with bounds

     runs in O(N log N)

     l; lower bound vector
     u; upper bound vector
     w; value vector
     N; size
     x; solution (if pointer different from NULL)

  */
  // map weights to value-index datastruct
  std::vector<associate, std::allocator<associate> > iw;
  NUM C = 0.0; // capacity
  NUM O = 0.0; // optimum
  iw.resize(N);
  for (unsigned int i = 0; i < N; ++i) {
    iw[i] = associate(w[i],i);
    C += l[i]; // compute lower bound sum
    O += l[i]*w[i]; // fill pack with minimum amounts
    if (x != NULL) x[i] = l[i];
  }
  // sort weights in non decreasing order
  sort(iw.begin(),iw.end());
    
  // fill variables up to capacity in reverse sorted order
  std::vector<associate, std::allocator<associate> >::iterator i = iw.end()-1;
  assert (C <= 1.0); // check feasibility
  while (C < 1.0) {
    NUM d = u[i->index]-l[i->index]; // remaining item's capacity 
    NUM r = 1.0 - C; // reamining capacity
    if (d < r) { // if there's still room to pack item completely, do it
      O += d*i->value;
      C += d;
      if (x != NULL) x[i->index] += d;
    } else { // otherwise, get as much as possible
      O += r*i->value;
      C += r;
      if (x != NULL) x[i->index] += r;
    }
    if (i == iw.begin()) break;
    --i;
  }
    
  return O;
}

NUM heapknapsack(NUM l[], NUM u[], NUM w[], unsigned int N, NUM* x) {
  /* Greedy implementation of the continuous knapsack problem with bounds
     with partial heap sorting

     runs in O(N log N)

     l; lower bound vector
     u; upper bound vector
     w; value vector
     N; size
     x; solution (if pointer different from NULL)

  */
  // map weights to value-index datastruct
  std::vector<associate, std::allocator<associate> > iw;
  NUM C = 0.0; // capacity
  NUM O = 0.0; // optimum
  iw.resize(N);
  for (unsigned int i = 0; i < N; ++i) {
    iw[i] = associate(w[i],i);
    C += l[i]; // compute lower bound sum
    O += l[i]*w[i]; // fill pack with minimum amounts
    if (x != NULL) x[i] = l[i];
  }
  
  // create a heap structure
  std::make_heap(iw.begin(),iw.end());
  
  // fill variables up to capacity by decreasing weight
  std::vector<associate, std::allocator<associate> >::iterator it = iw.end();
  assert (C <= 1.0); // check feasibility
  while (C < 1.0) {
    std::pop_heap(iw.begin(),it); // move heaviest element to the end of the vector
    --it;
    NUM d = u[it->index]-l[it->index]; // remaining item's capacity 
    NUM r = 1.0 - C; // reamining capacity
    if (d < r) { // if there's still room to pack item completely, do it
      O += d*it->value;
      C += d;
      if (x != NULL) x[it->index] += d;
    } else { // otherwise, get as much as possible
      O += r*it->value;
      C += r;
      if (x != NULL) x[it->index] += r;
    }
    if (it == iw.begin()) break;
  }
    
  return O;
}


NUM viterbi(NUM* A, NUM* B, NUM* P, unsigned int* O, unsigned int N, unsigned int M, unsigned int K, unsigned int* qstar, int debug_mode, int opt) {
  /* Computes Most Likely State Sequence (homogeneous models)
	 
     A; transition probability matrix
     B; emission probability matrix
     P; prior probability matrix
     O; observation sequence
     N; number of tates
     M; number of outcomes
     K; number of observations
     qstar; most likely state sequence
     debug_mode; output extra information
     opt; optimization algorithm
  */
    
  unsigned int qbest;
  NUM *phi = new NUM[K*N]; // decision matrix */
  NUM *d = new NUM[N];
  NUM *dt = new NUM[N];
    
  /* initial values */
  for (unsigned int i=0; i < N; ++i)
    d[i] = 0.0;
    
  NUM *l; // local transition credal set (Al=Au)
  // For t=K to t=2
  for (unsigned int t=K-1; t>0; --t) {
    for (unsigned int i=0; i < N; ++i) {
      dt[i] = -10000.0;
      l = A + i*N; // select local transition credal set A_i
      /* d_t(i) = max_q
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */
      for (unsigned int q=0; q < N; ++q) {
	/* compute
	   A[i][q]*B[q][O[t]]d_{t+1}(q) */
	NUM m = log(l[q]) + d[q] + log(B[O[t]*N+q]);
                
	if (m > dt[i]) {
	  dt[i] = m;
	  phi[t+i*K] = q;         
	}
      }
            
    }
    // update values
    for (unsigned int i=0; i < N; ++i)
      d[i] = dt[i];
        
  }
    
  // final update
  /* Q = max
     p[j]*B[j][O[0]] d_2(j) */
  qbest = 0;
  NUM Q = -10000.0;
  for (unsigned int q=0; q < N; ++q) {
    /*  compute 
	P[q]*B[q][O[t]]d_2(q) */
    NUM m = log(P[q]) + d[q] + log(B[q+N*O[0]]);
        
    if (m > Q) {
      Q = m;
      qbest = q;
    }
  }                           
    
  // MPE backtracking    
  qstar[0] = qbest;
    
  for (unsigned int t=1; t<K; ++t) {
    qbest = phi[t+qbest*K];
    qstar[t] = qbest;
  }
    
  // free memory
  delete [] d;
  delete [] dt;
  delete [] phi;    
	
  return Q;
}


std::pair<NUM,NUM> likelihood(NUM* Al, NUM* Au, NUM* Bl, NUM* Bu, NUM* Pl, NUM* Pu, unsigned int* O, unsigned int N, unsigned int M, unsigned int K, int debug_mode, int opt) {
  /* Log-Likelihood Computation (homogenous models)
	 
     Al; lower transition probability matrix
     Au; upper transition probability matrix
     Bl; lower emission probability matrix
     Bu; upper emission probability matrix
     Pl; lower prior probability matrix
     Pu; upper emission probability matrix
     O;  observation sequence
     N; number of tates
     M; number of outcomes
     K; number of observations
     debug_mode; output extra information
     opt; optimization algorithm
  */
	
  NUM *b = new NUM[N];
  NUM *w = new NUM[N];
  std::pair<NUM,NUM> ll(0.0, 0.0); // likelihood interval

	
  // Compute min loglikelihood
  /* initial values */
  for (unsigned int i=0; i < N; ++i)
    b[i] = 1.0;
	
  NUM *l, *u; // local transition credal set
  // For t=NUM to t=2
  for (unsigned int t=K-1; t>0; --t) {
    for (unsigned int j=0; j < N; ++j)
      w[j] = -b[j]*Bl[O[t]*N+j];
    NUM c = 0.0; // normalization scaling factor c_t
    for (unsigned int i=0; i < N; ++i) {
      l = Al + i*N; // select local transition credal set A_i
      u = Au + i*N;
      /* b_t(i) = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */
      b[i] = -knapsack(l,u,w,N,opt);
      c += b[i];
    }
    ll.first += log(c);
    for (unsigned int i=0; i < N; ++i)
      b[i] /= c;
  }
  // final update
  /* ll = max \sum_{j=1}^N
     p[j]*B[j][O[0]] b_2(j) */
  for (unsigned int j=0; j < N; ++j)
    w[j] = -b[j]*Bl[O[0]*N+j];
  /* b_0 = max_p \sum_{j=1}^N
     A[i][j]*B[j][O[t]]b_{t+1}(j) */
  ll.first += log(-knapsack(Pl,Pu,w,N,opt));
	
	
  // Compute max loglikelihood
  /* initial values */
  for (unsigned int i=0; i < N; ++i)
    b[i] = 1.0;
	
  // For t=K to t=2
  for (unsigned int t=K-1; t>0; --t) {
    for (unsigned int j=0; j < N; ++j)
      w[j] = b[j]*Bu[O[t]*N+j];
    NUM c = 0.0; // normalization scaling factor c_t
    for (unsigned int i=0; i < N; ++i) {
      l = Al + i*N; // select local transition credal set A_i
      u = Au + i*N;
      /* b_t(i) = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */
      b[i] = knapsack(l,u,w,N,opt);
      c += b[i];
    }
    ll.second += log(c);
    for (unsigned int i=0; i < N; ++i)
      b[i] /= c;
  }
  // final update
  /* ll = max \sum_{j=1}^N
     p[j]*B[j][O[0]] b_2(j) */
  for (unsigned int j=0; j < N; ++j)
    w[j] = b[j]*Bu[O[0]*N+j];
	
  /* b_0 = max_p \sum_{j=1}^N
     A[i][j]*B[j][O[t]]b_{t+1}(j) */
  ll.second += log(knapsack(Pl,Pu,w,N,opt));

  // free memory
  delete [] b;
  delete [] w;
    
  return ll;
}

std::pair<NUM,NUM> likelihood2(NUM** Al, NUM** Au, NUM** Bl, NUM** Bu, NUM* Pl, NUM* Pu, unsigned int* O, unsigned int* N, unsigned int* M, unsigned int K, int debug_mode, int opt) {
  /* Log-Likelihood Computation (non-homogenous models)
	 
     Al; lower transition probability matrix
     Au; upper transition probability matrix
     Bl; lower emission probability matrix
     Bu; upper emission probability matrix
     Pl; lower prior probability matrix
     Pu; upper emission probability matrix
     O;  observation sequence
     N; number of states
     M; number of outcomes
     K; number of observations
     debug_mode; output extra information
     opt; optimization algorithm
  */
	
  NUM** b = new NUM*[K]; // belief variables
  std::pair<NUM,NUM> ll(0.0, 0.0); // likelihood interval

  /* initialize b's */
  for (unsigned int t=0; t < K; ++t) b[t] = new NUM[N[t]];

  // Compute min loglikelihood
  /* initial values */
  for (unsigned int i=0; i < N[K-1]; ++i) b[K-1][i] = 1.0;
	
  NUM *l, *u; // local transition credal set
  NUM *w; // local weights
  // For t=K to t=2
  for (unsigned int t=K-1; t>0; --t) {
    w = new NUM[N[t]]; // weights w_t(s_t) = b_t[j][O[t]]*b_{t+1}(j)
    for (unsigned int j=0; j < N[t]; ++j)
      w[j] = -b[t][j]*Bl[t][O[t]+j*M[t]];
    NUM c = 0.0; // normalization scaling factor c_t
    for (unsigned int i=0; i < N[t-1]; ++i) {
      l = Al[t-1] + i*N[t]; // select local transition credal set A_i
      u = Au[t-1] + i*N[t];
      /* b_t(i) = min_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */
      b[t-1][i] = -knapsack(l,u,w,N[t],opt);
      c += b[t-1][i];
    }

    ll.first += log(c);
    for (unsigned int i=0; i < N[t-1]; ++i)
      b[t-1][i] /= c;
    delete [] w;// it'd be faster to allocate the maximum needed memory first and only deallocate in the end
  }
  // final update
  /* ll = min \sum_{j=1}^N  p[j]*B[j][O[0]] b_2(j) */
  w = new NUM[N[0]]; // weights w_1(s_1) = B_1[j][O[1]]*b_2(j)
  for (unsigned int j=0; j < N[0]; ++j)
    w[j] = -b[0][j]*Bl[0][O[0]+j*M[0]];
  /* b_0 = max_p \sum_{j=1}^N
     A[i][j]*B[j][O[t]]b_{t+1}(j) */
  ll.first += log(-knapsack(Pl,Pu,w,N[0],opt));
  delete [] w;
	
  // Compute max loglikelihood
  /* initial values */
  for (unsigned int i=0; i < N[K-1]; ++i) b[K-1][i] = 1.0;
	
  // For t=K to t=2
  for (unsigned int t=K-1; t>0; --t) {
    w = new NUM[N[t]]; // weights w_t(s_t) = b_t[j][O[t]]*b_{t+1}(j)
    for (unsigned int j=0; j < N[t]; ++j)
      w[j] = b[t][j]*Bu[t][O[t]+j*M[t]];
    NUM c = 0.0; // normalization scaling factor c_t
    for (unsigned int i=0; i < N[t-1]; ++i) {
      l = Al[t-1] + i*N[t]; // select local transition credal set A_i
      u = Au[t-1] + i*N[t];
      /* b_t(i) = min_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */
      b[t-1][i] = knapsack(l,u,w,N[t],opt);
      c += b[t-1][i];
    }

    ll.second += log(c);
    for (unsigned int i=0; i < N[t-1]; ++i)
      b[t-1][i] /= c;
    delete [] w;// it'd be faster to allocate the maximum needed memory first and only deallocate in the end
  }
  // final update
  /* ll = min \sum_{j=1}^N  p[j]*B[j][O[0]] b_2(j) */
  w = new NUM[N[0]]; // weights w_1(s_1) = B_1[j][O[1]]*b_2(j)
  for (unsigned int j=0; j < N[0]; ++j)
    w[j] = b[0][j]*Bu[0][O[0]+j*M[0]];
  /* b_0 = max_p \sum_{j=1}^N
     A[i][j]*B[j][O[t]]b_{t+1}(j) */
  ll.second += log(knapsack(Pl,Pu,w,N[0],opt));

  // free memory
  for (unsigned int t=0; t < K; ++t) delete[] b[t];
  delete [] b;
  delete [] w;
    
  return ll;
}


NUM viterbi2(NUM** A, NUM** B, NUM* P, unsigned int* O, unsigned int* N, unsigned int* M, unsigned int K, unsigned int* qstar, int debug_mode, int opt) {
  /* Computes Most Likely State Sequence (non-homogenous models)
     Fills in array qstart containing an optimal solution

     A; transition probability matrices
     B; emission probability matrices
     P; prior probability matrix
     O;  observation sequence
     N; array of number of states
     M; array of number of outcomes
     K; number of observations
     qstar; most likely state sequence (output)
     debug_mode; output extra information
     opt; optimization algorithm
  */
    
  unsigned int qbest;
  unsigned int** phi = new unsigned int*[K]; // decision matrix 
  NUM** d = new NUM*[K]; // propagation variables 
    
  /* initial values */
  for (unsigned int k=0; k < K; ++k) {
    phi[k] = new unsigned int[N[k]];
    d[k] = new NUM[N[k]];
    for (unsigned int i=0; i < N[k]; ++i)
      d[k][i] = 0.0;
  }
    
  NUM *l; //  transition distribution
  // For t=K to t=2
  for (unsigned int t=K-1; t > 0; --t) {
    for (unsigned int i=0; i < N[t-1]; ++i) {
      d[t-1][i] = -1e+10; // -infitinity 
      l = A[t-1] + i*N[t]; // select local transition credal set A_i
      /* d_{t-1}(i) = max_{q \in N_t} A_t[i][q]*B_t[q][O[t]]*d_t(q) */
      for (unsigned int q=0; q < N[t]; ++q) {
	/* compute A_t[i][q]*B_t[q][O[t]]*d_t(q) */
	NUM m = log(l[q]) + d[t][q] + log(B[t][O[t]+q*M[t]]);
                
	if (m > d[t-1][i]) {
	  d[t-1][i] = m;
	  phi[t][i] = q;         
	}
      }
            
    }        
  }
    
  // final update
  /* Q = max_{j \in N_0} p[j]*B[j][O[0]]*d_0(j) */
  qbest = 0;
  NUM Q = -10000.0;
  for (unsigned int q=0; q < N[0]; ++q) {
    /* compute P[q]*B[q][O[t]]*d_0(q) */
    NUM m = log(P[q]) + d[0][q] + log(B[0][O[0]+q*M[0]]);
        
    if (m > Q) {
      Q = m;
      qbest = q;
    }
  }                           
    
  // MPE backtracking    
  qstar[0] = qbest;
    
  for (unsigned int t=1; t<K; ++t) {
    qbest = phi[t][qbest];
    qstar[t] = qbest;
  }
    
  // free memory
  for (unsigned int k=0; k < K; ++k) { 
      delete [] d[k];
      delete [] phi[k];    
  }
  delete [] d;
  delete [] phi;   
	
  return Q;
}

std::pair<NUM,NUM> filtering(unsigned q, NUM* Al, NUM* Au, NUM* Bl, NUM* Bu, NUM* Pl, NUM* Pu, unsigned int* O, unsigned int N, unsigned int M, unsigned int K, double eps, int debug_mode, int opt) {
  /* Filtering (homogenous models)
	 
     q;  state of current hidden variable to predict
     Al; lower transition probability matrix
     Au; upper transition probability matrix
     Bl; lower emission probability matrix
     Bu; upper emission probability matrix
     Pl; lower prior probability matrix
     Pu; upper emission probability matrix
     O;  observation sequence
     N; number of tates
     M; number of outcomes
     K; number of observations
     eps; binary search precision
     debug_mode; output extra information
     opt; optimization algorithm
  */
	
  NUM *g = new NUM[N];
  NUM *w = new NUM[N];
  std::pair<NUM,NUM> prob(0.0, 0.0); // likelihood interval

  NUM *l, *u; // local transition credal set
  
  NUM left = 0.0, right = 1.0;
  NUM k = 0.0;

  // lower probability

  if (debug_mode) { std::cout << std::endl; }  
  
  while (right - left > eps)
    {
      k = (left+right)/2;

      /* initial values */
      for (unsigned int i=0; i < N; ++i)
	if (i==q)
	  g[i] = 1.0-k;
	else
	  g[i] = -k;
	
      // For t=K to t=2
      for (unsigned int t=K-1; t>0; --t) {
	for (unsigned int j=0; j < N; ++j)
	  {
	    if (g[j] > 0)
	      w[j] = -g[j]*Bl[O[t]*N+j];
	    else
	      w[j] = -g[j]*Bu[O[t]*N+j];
	  }
	for (unsigned int i=0; i < N; ++i) {
	  l = Al + i*N; // select local transition credal set A_i
	  u = Au + i*N;
	  /* b_t(i) = max_p \sum_{j=1}^N
	     A[i][j]*B[j][O[t]]b_{t+1}(j) */
	  g[i] = -knapsack(l,u,w,N,opt);
	}
      }
      // final update
      /* ll = max \sum_{j=1}^N
	 p[j]*B[j][O[0]] b_2(j) */
      for (unsigned int j=0; j < N; ++j)
	{
	  if (g[j] > 0)
	    w[j] = -g[j]*Bl[O[0]*N+j];
	  else
	    w[j] = -g[j]*Bu[O[0]*N+j];
	}
      /* b_0 = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */  
      prob.first = -knapsack(Pl,Pu,w,N,opt);

      if (debug_mode) {
	std::cout << std::fixed << " [" << left << ", " << right << "] ";// << std::endl;
	std::cout << prob.first << "\n";
      }

  
      // check probability signal
      if (prob.first > 0.0)
	left = k; // update interval left bound
      else
	right = k; // update interval right bound
    }

  prob.first = (left+right)/2.0;
  
  if (debug_mode) { std::cout << std::endl; }


  // upper probability

  left = 0.0; right = 1.0;
  
  while (right - left > eps)
    {
      k = (left+right)/2;

      /* initial values */
      for (unsigned int i=0; i < N; ++i)
	if (i==q)
	  g[i] = 1.0-k;
	else
	  g[i] = -k;
	
      // For t=K to t=2
      for (unsigned int t=K-1; t>0; --t) {
	for (unsigned int j=0; j < N; ++j)
	  {
	    if (g[j] > 0)
	      w[j] = g[j]*Bu[O[t]*N+j];
	    else
	      w[j] = g[j]*Bl[O[t]*N+j];
	  }
	for (unsigned int i=0; i < N; ++i) {
	  l = Al + i*N; // select local transition credal set A_i
	  u = Au + i*N;
	  /* b_t(i) = max_p \sum_{j=1}^N
	     A[i][j]*B[j][O[t]]b_{t+1}(j) */
	  g[i] = knapsack(l,u,w,N,opt);
	}
      }
      // final update
      /* ll = max \sum_{j=1}^N
	 p[j]*B[j][O[0]] b_2(j) */
      for (unsigned int j=0; j < N; ++j)
	{
	  if (g[j] > 0)
	    w[j] = g[j]*Bu[O[0]*N+j];
	  else
	    w[j] = g[j]*Bl[O[0]*N+j];
	}
      /* b_0 = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */  
      prob.second = knapsack(Pl,Pu,w,N,opt);

      if (debug_mode) {
	std::cout << std::fixed << " [" << left << ", " << right << "] ";// << std::endl;
	std::cout << prob.second << "\n";
      }

  
      // check probability signal
      if (prob.second > 0.0)
	left = k; // update interval left bound
      else
	right = k; // update interval right bound
    }

  prob.second = (left+right)/2.0;
  
  if (debug_mode) { std::cout << std::endl; }  
  
  // free memory
  delete [] g;
  delete [] w;
  return prob;

}

std::pair<NUM,NUM> filtering2(unsigned q, NUM** Al, NUM** Au, NUM** Bl, NUM** Bu, NUM* Pl, NUM* Pu, unsigned int* O, unsigned int* N, unsigned int* M, unsigned int K, double eps, int debug_mode, int opt) {
  /* Filtering (homogenous models)
	 
     q;  state of current hidden variable to predict
     Al; lower transition probability matrix
     Au; upper transition probability matrix
     Bl; lower emission probability matrix
     Bu; upper emission probability matrix
     Pl; lower prior probability matrix
     Pu; upper emission probability matrix
     O;  observation sequence
     N; number of states
     M; number of outcomes
     K; number of observations
     eps; binary search precision
     debug_mode; output extra information
     opt; optimization algorithm
  */
	
  NUM** g = new NUM*[K];
  /* initialize g */
  for (unsigned int t=0; t < K; ++t) g[t] = new NUM[N[t]];

  NUM *w;
  std::pair<NUM,NUM> prob(0.0, 0.0); // likelihood interval

  NUM *l, *u; // local transition credal sets
  
  NUM left = 0.0, right = 1.0;
  NUM k = 0.0;

  // lower probability

  if (debug_mode) { std::cout << std::endl; }  
  
  while (right - left > eps)
    {
      k = (left+right)/2;

      /* initial values */
      for (unsigned int i=0; i < N[K-1]; ++i)
	if (i==q)
	  g[K-1][i] = 1.0-k;
	else
	  g[K-1][i] = -k;
	
      // For t=K to t=2
      for (unsigned int t=K-1; t>0; --t) {
	w = new NUM[N[t]]; // weights 
	for (unsigned int j=0; j < N[t]; ++j)
	  {
	    if (g[t][j] > 0)
	      w[j] = -g[t][j]*Bl[t][O[t]+j*M[t]];
	    else
	      w[j] = -g[t][j]*Bu[t][O[t]+j*M[t]];
	  }
	for (unsigned int i=0; i < N[t]; ++i) {
	  l = Al[t-1] + i*N[t]; // select local transition credal set A_i
	  u = Au[t-1] + i*N[t];
	  /* b_t(i) = max_p \sum_{j=1}^N
	     A[i][j]*B[j][O[t]]b_{t+1}(j) */
	  g[t-1][i] = -knapsack(l,u,w,N[t],opt);
	}
	delete [] w;// it'd be faster to allocate the maximum needed memory first and only deallocate in the end
      }
      // final update
      /* ll = max \sum_{j=1}^N
	 p[j]*B[j][O[0]] b_2(j) */
      w = new NUM[N[0]]; // weights 
      for (unsigned int j=0; j < N[0]; ++j)
	{
	  if (g[0][j] > 0)
	    w[j] = -g[0][j]*Bl[0][O[0]+j*M[0]];
	  else
	    w[j] = -g[0][j]*Bu[0][O[0]+j*M[0]];
	}
      /* b_0 = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */  
      prob.first = -knapsack(Pl,Pu,w,N[0],opt);
      delete [] w;

      if (debug_mode) {
	std::cout << std::fixed << " [" << left << ", " << right << "] ";// << std::endl;
	std::cout << prob.first << "\n";
      }

  
      // check probability signal
      if (prob.first > 0.0)
	left = k; // update interval left bound
      else
	right = k; // update interval right bound
    }

  prob.first = (left+right)/2.0;
  
  if (debug_mode) { std::cout << std::endl; }


  // upper probability

  left = 0.0; right = 1.0;
  
  if (debug_mode) { std::cout << std::endl; }  
  
  while (right - left > eps)
    {
      k = (left+right)/2;

      /* initial values */
      for (unsigned int i=0; i < N[K-1]; ++i)
	if (i==q)
	  g[K-1][i] = 1.0-k;
	else
	  g[K-1][i] = -k;
	
      // For t=K to t=2
      for (unsigned int t=K-1; t>0; --t) {
	w = new NUM[N[t]]; // weights 
	for (unsigned int j=0; j < N[t]; ++j)
	  {
	    if (g[t][j] > 0)
	      w[j] = g[t][j]*Bu[t][O[t]+j*M[t]];
	    else
	      w[j] = g[t][j]*Bl[t][O[t]+j*M[t]];
	  }
	for (unsigned int i=0; i < N[t]; ++i) {
	  l = Al[t-1] + i*N[t]; // select local transition credal set A_i
	  u = Au[t-1] + i*N[t];
	  /* b_t(i) = max_p \sum_{j=1}^N
	     A[i][j]*B[j][O[t]]b_{t+1}(j) */
	  g[t-1][i] = knapsack(l,u,w,N[t],opt);
	}
	delete [] w;// it'd be faster to allocate the maximum needed memory first and only deallocate in the end
      }
      // final update
      /* ll = max \sum_{j=1}^N
	 p[j]*B[j][O[0]] b_2(j) */
      w = new NUM[N[0]]; // weights 
      for (unsigned int j=0; j < N[0]; ++j)
	{
	  if (g[0][j] > 0)
	    w[j] = g[0][j]*Bu[0][O[0]+j*M[0]];
	  else
	    w[j] = g[0][j]*Bl[0][O[0]+j*M[0]];
	}
      /* b_0 = max_p \sum_{j=1}^N
	 A[i][j]*B[j][O[t]]b_{t+1}(j) */  
      prob.second = knapsack(Pl,Pu,w,N[0],opt);
      delete [] w;

      if (debug_mode) {
	std::cout << std::fixed << " [" << left << ", " << right << "] ";// << std::endl;
	std::cout << prob.first << "\n";
      }

  
      // check probability signal
      if (prob.second > 0.0)
	left = k; // update interval left bound
      else
	right = k; // update interval right bound
    }

  prob.second = (left+right)/2.0;
  
  if (debug_mode) { std::cout << std::endl; }
  
  // free memory
  for (unsigned int t=0; t < K; ++t) delete[] g[t];
  delete [] g;
  return prob;

}
