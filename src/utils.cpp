/* 
 * imprecise hidden Markov model
 *
 * Utility functions
 * 
 * See COPYRIGHT and VERSION files for copyright and version info
 *
 */

#include "utils.h"

// input

bool load_model2(NUM **&Al, NUM **&Au, NUM **&Bl, NUM **&Bu, NUM *&Pl, NUM *&Pu, unsigned int* &N, unsigned int* &M, unsigned int &K, char* filename) {
  // Load non homogeneous model from file
  // return 0 on success, 1 on failure

  std::ifstream ifs(filename); // open file stream
	
  std::string line;

  // read model datafile
  if (ifs.is_open())
    {
      unsigned state = 0; // initial state
      unsigned sizeA = 0; // size of matrix A
      unsigned sizeB = 0; // size of matrix B
      NUM value;
      while(std::getline(ifs, line)) {
	if ((line.size() > 0) && !(line[0] == '#')) { 
	  // ignore lines starting with #
	  switch(state) {
	    
	  case 0:
	    { // read K
	      std::istringstream iss(line, std::istringstream::in);
	      iss >> K;
	      ++state;
	    }
	    break;
	    
	  case 1:
	    { // read N
	      N = new unsigned int[K];
	      if (N == NULL) {
		std::cerr << "Could not initialize number of states array!" << std::endl;
		return 1; // FAILURE
	      } 
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int v, n=0;
	      while( iss >> v ) {
		if (n >= K) {
		  std::cerr << "Format mismatch: Dimension mismatch in number of states array!" << std::endl;
		  return 1;	// FAILURE	
		}
		N[n] = v;
		if (n) sizeA += N[n-1]*v; // num of parameters P(s_n|s_{n-1})
		n++;
	      } 
	      ++state;
	    }
	    break;

	  case 2:
	    { // read M
	      M = new unsigned int[K];
	      if (M == NULL) {
		std::cerr << "Could not initialize number of symbols array!" << std::endl;
		return 1; // FAILURE
	      } 
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int v, n=0;
	      while( iss >> v ) {
		if (n >= K) {
		  std::cerr << "Format mismatch: Dimension mismatch in number of symbols array!" << std::endl;
		  return 1;	// FAILURE	
		}
		M[n] = v;
		sizeB += v*N[n]; // num of parameters P(o_n|s_n)
		n++;
	      }
	      ++state;
	    }
	    break;

	  case 3:
	    {
	      // initialize matrices
	      Al = new NUM*[K-1];
	      Au = new NUM*[K-1];
	      if (Al == NULL || Au == NULL) {
		std::cerr << "Could not allocate transition matrices!" << std::endl;
		return 1; // FAILURE
	      } 

	     for (unsigned n=1; n<K; n++) {
		// initialize matrices
		Al[n-1] = new NUM[N[n]*N[n-1]];
		Au[n-1] = new NUM[N[n]*N[n-1]];
		if (Al[n-1] == NULL || Au[n-1] == NULL) {
		  std::cerr << "Could not allocate transition matrices beyond step " << n << "!" << std::endl;
		  return 1; // FAILURE
		}
	     }
		           
	      // read lower transition probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0, k = 1, t = 0;
	      while( iss >> value ) {
		/* read lower transition matrix */
		if (t >= sizeA) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower transition matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Al[k-1][n] = value;
		n++;
		if (n >= N[k]*N[k-1]) {
		  n = 0; k++;
		}
		t++;
	      }
	      if (t < sizeA) {
		std::cerr << "Format mismatch: Dimension mismatch in lower transition matrix!" << std::endl;
		return 1; // FAILURE		
	      }						      
	      ++state;
	    }
	    break;

	  case 4:
	    {
	      // read upper transition probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0, k = 1, t = 0;
	      while( iss >> value ) {
		/* read lower transition matrix */
		if (t >= sizeA) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper transition matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Au[k-1][n] = value;
		n++;
		if (n >= N[k]*N[k-1]) {
		  n = 0; k++;
		}
		t++;
	      }
	      if (t < sizeA) {
		std::cerr << "Format mismatch: Dimension mismatch in upper transition matrix!" << std::endl;
		return 1; // FAILURE		
	      }						      
	      ++state;
	    }
	    break;
	      
	  case 5:
	    {
	      // initialize matrices
	      Bl = new NUM*[K];
	      Bu = new NUM*[K];
	      if (Bl == NULL || Bu == NULL) {
		std::cerr << "Could not allocate emission matrices!" << std::endl;
		return 1; // FAILURE
	      } 

	     for (unsigned n=0; n<K; n++) {
		// initialize matrices
		Bl[n] = new NUM[M[n]*N[n]];
		Bu[n] = new NUM[M[n]*N[n]];
		if (Bl[n] == NULL || Bu[n] == NULL) {
		  std::cerr << "Could not allocate emission matrices beyond step " << n << "!" << std::endl;
		  return 1; // FAILURE
		}
	     }

	      // read lower emission probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0, t = 0, k = 0;
	      while( iss >> value ) {
		/* read lower emission matrix */
		if (t >= sizeB) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower emission matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Bl[k][n] = value;
		n++; 
		if (n >= M[k]*N[k]) {
		  n = 0; k++;
		}
		t++;
	      }
	      if (t < sizeB) {
		std::cerr << "Format mismatch: Dimension mismatch in lower emission matrix!" << std::endl;
		return 1; // FAILURE		
	      }		
	      ++state;
	    }
	    break;

	  case 6:
	    {
	      // read upper emission probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0, t = 0, k = 0;
	      while( iss >> value ) {
		/* read lower emission matrix */
		if (t >= sizeB) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper emission matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Bu[k][n] = value;
		n++; 
		if (n >= M[k]*N[k]) {
		  n = 0; k++;
		}
		t++;
	      }
	      if (t < sizeB) {
		std::cerr << "Format mismatch: Dimension mismatch in upper emission matrix!" << std::endl;
		return 1; // FAILURE		
	      }		
	      ++state;
	    }
	    break;

	  case 7:
	    {
	      // initialize matrices
	      Pl = new NUM[N[0]];
	      Pu = new NUM[N[0]];        
	      if (Pl == NULL || Pu == NULL) {
		std::cerr << "Could not initialize prior matrices!" << std::endl;
		return 1; // FAILURE	
	      }
	      // read lower prior probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read lower prior matrix */
		if (n >= N[0]) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower prior matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Pl[n] = value;
		n++;
	      }
	      if (n < N[0]) {
		std::cerr << "Format mismatch: Dimension mismatch in lower prior matrix!" << std::endl;
		return 1; // FAILURE		
	      }					
	      ++state;
	    }
	    break;

	  case 8:
	    {	
	      // read upper prior probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read upper prior matrix */
		if (n >= N[0]) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper prior matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Pu[n] = value;
		n++;
	      }
	      if (n < N[0]) {
		std::cerr << "Format mismatch: Dimension mismatch in upper prior matrix!" << std::endl;
		return 1;  // FAILURE	
	      }					
	      ++state;
	    }
	    break;

	  } // end of automata
	}
      } // no more lines to read
	/* check if model was correctly loaded: was each state reached? */
      if (state < 1)
	std::cerr << "Format mismatch: Expected integer: number of steps!" << std::endl;

      if (state < 2)
	std::cerr << "Format mismatch: Expected integer: number of states!" << std::endl;

      if (state < 3)
	std::cerr << "Format mismatch: Expected integer: Number of outcomes!" << std::endl;

      if (state < 4)
	std::cerr << "Format mismatch: Lower transition matrix not found!" << std::endl;

      if (state < 5)
	std::cerr << "Format mismatch: Upper transition matrix not found!" << std::endl;

      if (state < 6)
	std::cerr << "Format mismatch: Lower emission matrix not found!" << std::endl;

      if (state < 7)
	std::cerr << "Format mismatch: Upper transition matrix not found!" << std::endl;

      if (state < 8)
	std::cerr << "Format mismatch: Lower prior matrix not found!" << std::endl;

      if (state < 9)
	std::cerr << "Format mismatch: Upper prior matrix not found!" << std::endl;

    } else { // file couldn't be opened
    std::cerr << "Error: file " << filename << " could not be opened" << std::endl;
    return 1; // FAILURE
  }

  /* close input file stream */
  ifs.close();

  return 0; // OK

}

bool load_model(NUM *&Al, NUM *&Au, NUM *&Bl, NUM *&Bu, NUM *&Pl, NUM *&Pu, unsigned int &N, unsigned int &M, char* filename) {
  // Load homogeneous model from file
  // return 0 on success, 1 on failure

  std::ifstream ifs(filename); // open file stream
	
  std::string line;

  // read model datafile
  if (ifs.is_open())
    {
      int state = 0; // initial state
      NUM value;
      while(std::getline(ifs, line)) {
	if ((line.size() > 0) && !(line[0] == '#')) { 
	  // ignore lines starting with #
	  switch(state) {
	    
	  case 0:
	    {
	      // read N
	      std::istringstream iss(line, std::istringstream::in);
	      iss >> N;
	      ++state;
	    }
	    break;

	  case 1:
	    {
	      // read M
	      std::istringstream iss(line, std::istringstream::in);
	      iss >> M;
	      ++state;
	    }
	    break;

	  case 2:
	    {
	      // initialize matrices
	      Al = new NUM[N*N];
	      Au = new NUM[N*N];
	      if (Al == NULL || Au == NULL) {
		std::cerr << "Could not initialize transition matrix!" << std::endl;
		return 1; // FAILURE
	      } 
	      Bl = new NUM[N*M];
	      Bu = new NUM[N*M];
	      if (Bl == NULL || Bu == NULL) {
		std::cerr << "Could not initialize emission matrix!" << std::endl;
		return 1; // FAILURE
	      }
	      Pl = new NUM[N];
	      Pu = new NUM[N];        
	      if (Pl == NULL || Pu == NULL) {
		std::cerr << "Could not initialize prior matrix!" << std::endl;
		return 1; // FAILURE	
	      }
        
	      // read lower transition probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read lower transition matrix */
		if (n >= N*N) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower transition matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Al[n] = value;
		n++;
	      }
	      if (n < N*N) {
		std::cerr << "Format mismatch: Dimension mismatch in lower transition matrix!" << std::endl;
		return 1; // FAILURE		
	      }						      
	      ++state;
	    }
	    break;

	  case 3:
	    {
	      // read upper transition probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read upper transition matrix */
		// store value
		if (n >= N*N) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper transition matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		Au[n] = value;
		n++;
	      }
	      if (n < N*N) {
		std::cerr << "Format mismatch: Dimension mismatch in upper transition matrix!" << std::endl;
		return 1; // FAILURE		
	      }						      
	      ++state;
	    }
	    break;
	      
	  case 4:
	    {
	      // read lower emission probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      unsigned int i = 0, j = 0;
	      while( iss >> value ) {
		/* read lower emission matrix */
		if (n >= N*M) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower emission matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Bl[i] = value;
		n++; 
		i+=N; // read matrix transpose
		if (i >= N*M) {
		  j++;
		  i = j;
		}
	      }
	      if (n < N*M) {
		std::cerr << "Format mismatch: Dimension mismatch in lower emission matrix!" << std::endl;
		return 1; // FAILURE		
	      }		
	      ++state;
	    }
	    break;

	  case 5:
	    {
	      // read upper emission probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      unsigned int i=0, j=0;
	      while( iss >> value ) {
		/* read upper emission matrix */
		if (n >= N*M) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper emission matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Bu[i] = value;
		n++; 
		i+=N; // read matrix transpose
		if (i >= N*M) {
		  j++;
		  i = j;
		}				
	      }
	      if (n < N*M) {
		std::cerr << "Format mismatch: Dimension mismatch in upper emission matrix!" << std::endl;
		return 1; // FAILURE		
	      }					
	      ++state;
	    }
	    break;

	  case 6:
	    {
	      // read lower prior probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read lower prior matrix */
		if (n >= N) {
		  std::cerr << "Format mismatch: Dimension mismatch in lower prior matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Pl[n] = value;
		n++;
	      }
	      if (n < N) {
		std::cerr << "Format mismatch: Dimension mismatch in lower prior matrix!" << std::endl;
		return 1; // FAILURE		
	      }					
	      ++state;
	    }
	    break;

	  case 7:
	    {	
	      // read upper prior probabilities
	      std::istringstream iss(line, std::istringstream::in);
	      unsigned int n = 0;
	      while( iss >> value ) {
		/* read upper prior matrix */
		if (n >= N) {
		  std::cerr << "Format mismatch: Dimension mismatch in upper prior matrix!" << std::endl;
		  return 1;	// FAILURE	
		}					
		// store value
		Pu[n] = value;
		n++;
	      }
	      if (n < N) {
		std::cerr << "Format mismatch: Dimension mismatch in upper prior matrix!" << std::endl;
		return 1;  // FAILURE	
	      }					
	      ++state;
	    }
	    break;

	  } // end of automata
	}
      } // no more lines to read
	/* check if model was correctly loaded: was each state reached? */
      if (state < 1)
	std::cerr << "Format mismatch: Expected integer: number of states!" << std::endl;

      if (state < 2)
	std::cerr << "Format mismatch: Expected integer: Number of outcomes!" << std::endl;

      if (state < 3)
	std::cerr << "Format mismatch: Lower transition matrix not found!" << std::endl;

      if (state < 4)
	std::cerr << "Format mismatch: Upper transition matrix not found!" << std::endl;

      if (state < 5)
	std::cerr << "Format mismatch: Lower emission matrix not found!" << std::endl;

      if (state < 6)
	std::cerr << "Format mismatch: Upper transition matrix not found!" << std::endl;

      if (state < 7)
	std::cerr << "Format mismatch: Lower prior matrix not found!" << std::endl;

      if (state < 8)
	std::cerr << "Format mismatch: Upper prior matrix not found!" << std::endl;

    } else { // file couldn't be opened
    std::cerr << "Error: file " << filename << " could not be opened" << std::endl;
    return 1; // FAILURE
  }

  /* close input file stream */
  ifs.close();

  return 0; // OK

}

bool load_data2(std::vector<sequence>& D, unsigned int K, unsigned int* M, char* filename) {
  // load observations from file (non-uniform model)
  // return 0 on success, 1 on failure
  std::ifstream dfs(filename); //  open file stream
  std::string line;    

  if (dfs.is_open())
    {
      unsigned int l=0;
      while (getline(dfs, line))
	{
	  l++;
	  if ((line.size() > 0) && !(line[0] == '#')) { // ignore lines starting with '#'
	    std::istringstream iss(line, std::istringstream::in);
            
	    unsigned int value;
	    unsigned int n=0;
	    sequence obs;
	    // read sequence
	    while( iss >> value )  {   
	      if (n >= K) { // sequence is too lengthy
		std::cerr << ">>> In line " << l << " of file " << filename << std::endl;
		std::cerr << "Invalid observation: given sequence contains more steps than model!" << std::endl;
		return 1;       // FAILURE
	      }

	      if (value >= M[n]) { // if value is greater than number of symbols
		std::cerr << ">>> In line " << l << " of file " << filename << std::endl;
		std::cerr << "Invalid observation: found value not in symbol range!" << std::endl;
		std::cerr << "Symbol found: " << value << std::endl;
		std::cerr << "Symbol range: " << 0 << "-" << (M[n]-1) << std::endl;
		return 1;	// FAILURE
	      }
	      obs.push_back(value);
	      n++;	      
	    }
	    if (obs.size() != K) // if sequence has more than 1 element, add it to dataset
	      std::cerr << "Invalid observation sequence: model and sequence length differ!" << std::endl;
	    else
	      D.push_back(obs);
	  }
	}

    } else { // file couldn't be opened
    std::cerr << "Error: file " << filename << " could not be opened" << std::endl;
    return 1;
  }	
	
  dfs.close(); // close data stream

  return 0;

}

bool load_data(std::vector<sequence>& D, unsigned int K, unsigned int M, char* filename) {
  // load observations from file
  // return 0 on success, 1 on failure
  std::ifstream dfs(filename); //  open file stream
  std::string line;    

  if (dfs.is_open())
    {
      unsigned int l=0;
      while (getline(dfs, line))
	{
	  l++;
	  if ((line.size() > 0) && !(line[0] == '#')) { // ignore lines starting with '#'
	    std::istringstream iss(line, std::istringstream::in);
            
	    unsigned int value;
	    sequence obs;
	    // read sequence
	    while( iss >> value )  {   
	      if (value >= M) { // if value is greater than emission matrix dimension
		std::cerr << ">>> In line " << l << " of file " << filename << std::endl;
		std::cerr << "Invalid observation: found observation value not in symbol range!" << std::endl;
		std::cerr << "Symbol found: " << value << std::endl;
		std::cerr << "Symbol range: " << 0 << "-" << (M-1) << std::endl;
		return 1;	// FAILURE
	      }
	      obs.push_back(value);
	    }
	    if (obs.size()) // if sequence has more than 1 element, add it to dataset
	      D.push_back(obs);
	  }
	}

    } else { // file couldn't be opened
    std::cerr << "Error: file " << filename << " could not be opened" << std::endl;
    return 1;
  }	
	
  dfs.close(); // close data stream

  return 0;

}

// output

void print_precise_probabilities (NUM **A, NUM **B, NUM *P, unsigned int* N, unsigned int* M, unsigned int K) 
{
    std::cout << "Prior probabilities" << std::endl;
    for (unsigned int i=0; i<N[0]; ++i)
      std::cout << "    pi(" << i << ")=" << P[i] << "    ";
    std::cout << std::endl;

    std::cout << "Transition probabilities" << std::endl;
    for (unsigned int n=1; n<K; ++n) {
      std::cout << "t=" << (n+1) << std::endl;
      for (unsigned int i=0; i<N[n-1]; ++i) {
	for (unsigned int j=0; j<N[n]; ++j) 
	  std::cout << "    a(" << j << "|" << i << ")=" << A[n-1][j+i*N[n]] << "    ";
	std::cout << std::endl;
      }	
    }
    std::cout << std::endl;

    std::cout << "Emission probabilities" << std::endl;
    for (unsigned int n=0; n<K; ++n) {
      std::cout << "t=" << (n+1) << std::endl;
      for (unsigned int i=0; i<N[n]; ++i) {
	for (unsigned int j=0; j<M[n]; ++j) 
	  std::cout << "    b(" << j << "|" << i << ")=" << B[n][j+i*M[n]] << "    ";
	std::cout << std::endl;
      }	
    }			
    std::cout << std::endl;
		
}

void print_probabilities2 (NUM **Al, NUM **Au, NUM **Bl, NUM **Bu, NUM *Pl, NUM *Pu, unsigned int* N, unsigned int* M, unsigned int K) 
{
    std::cout << "Prior probabilities" << std::endl;
    for (unsigned int i=0; i<N[0]; ++i)
      std::cout << "    pi(" << i << ") in [" << Pl[i] << "," << Pu[i] << "]    ";
    std::cout << (is_reachable(Pl,Pu,N[0])?"":"*** non reachable ***") << std::endl;    
    std::cout << std::endl;

    std::cout << "Transition probabilities" << std::endl;
    NUM *l, *u;
    for (unsigned int n=1; n<K; ++n) {
      std::cout << "t=" << (n+1) << std::endl;
      for (unsigned int i=0; i<N[n-1]; ++i) {
	for (unsigned int j=0; j<N[n]; ++j) 
	  std::cout << "    a(" << j << "|" << i << ") in [" << Al[n-1][j+i*N[n]] << "," << Au[n-1][j+i*N[n]] << "]    ";
	l = Al[n-1] + i*N[n];
	u = Au[n-1] + i*N[n];
	std::cout << (is_reachable(l,u,N[n])?"":"*** non reachable ***") << std::endl;    
	std::cout << std::endl;
      }	
    }
    std::cout << std::endl;

    std::cout << "Emission probabilities" << std::endl;
    for (unsigned int n=0; n<K; ++n) {
      std::cout << "t=" << (n+1) << std::endl;
      for (unsigned int i=0; i<N[n]; ++i) {
	for (unsigned int j=0; j<M[n]; ++j) 
	  std::cout << "    b(" << j << "|" << i << ") in [" << Bl[n][j+i*M[n]] << "," << Bu[n][j+i*M[n]] << "]    ";
	l = Bl[n] + i*M[n];
	u = Bu[n] + i*M[n];
	std::cout << (is_reachable(l,u,M[n])?"":"*** non reachable ***") << std::endl;    
	std::cout << std::endl;
      }	
    }			
    std::cout << std::endl;
		
}

void print_probabilities (NUM *Al, NUM *Au, NUM *Bl, NUM *Bu, NUM *Pl, NUM *Pu, unsigned int N, unsigned int M) 
{
    std::cout << "Prior probabilities" << std::endl;
    for (unsigned int i=0; i<N; ++i)
      std::cout << "    pi(" << i << ") in [" << Pl[i] << "," << Pu[i] << "]    ";
    std::cout << (is_reachable(Pl,Pu,N)?"":"*** non reachable ***") << std::endl;    
    std::cout << std::endl;
  
    std::cout << "Transition probabilities" << std::endl;
    NUM *l, *u;
    for (unsigned int i=0; i<N; ++i) {
      for (unsigned int j=0; j<N; ++j) 
	std::cout << "    a(" << j << "|" << i << ") in [" << Al[j+i*N] << "," << Au[j+i*N] << "]    ";
      l = Al + i*N;
      u = Au + i*N;
      std::cout << (is_reachable(l,u,N)?"":"*** non reachable ***") << std::endl;    
      std::cout << std::endl;
    }	
    std::cout << std::endl;

    std::cout << "Emission probabilities" << std::endl;
    for (unsigned int i=0; i<N; ++i) {
      NUM* l = new NUM[M];
      for (unsigned j=0; j<M; j++)
	l[j] = Bl[j*N+i]; 
      NUM* u = new NUM[M];
      for (unsigned j=0; j<M; j++)
	u[j] = Bu[j*N+i];       
      for (unsigned int j=0; j<M; ++j) 
	std::cout << "    b(" << j << "|" << i << ") in [" << l[j] << "," << u[j] << "]    ";
      std::cout << (is_reachable(l,u,M)?"":"*** non reachable ***") << std::endl;    
      std::cout << std::endl;
    }	
    std::cout << std::endl;    
}


void print_observations (std::vector<sequence> D)
{
  std::cout << "Observations" << std::endl;
  for (unsigned int s=0; s<D.size(); ++s) {
    std::cout << s << ". ";
    for (sequence::iterator o = D[s].begin(); o != D[s].end(); ++o)
      std::cout << *o << " ";
    std::cout << std::endl;
  }
}

// validation

bool is_reachable(NUM* l, NUM* u, NUM N) {
  // check whether a given hypercube set is reachable
  // input: arrays of lower and upper bounds l and u, dimension N
  for (unsigned i=0; i<N; i++) {
    // u_i + \sum_{j \neq i} l_j <= 1 <= l_i + \sum_{j \neq i} u_j
    NUM suml = (NUM)0.0, sumu = (NUM)0.0;
    for (unsigned j=0; j<N; j++) if (j != i) { suml += l[j]; sumu += u[j]; }
    if ((u[i] + suml > 1.0+1e-10) || (l[i] + sumu < 1.0-1e-10)) return false;
  }
  return true;
}

// time measurement

double getcputime(void) {
  // get total cpu time (in seconds) consumed by process
  struct timeval tim;        
  struct rusage ru;        
  getrusage(RUSAGE_SELF, &ru);        
  tim=ru.ru_utime;        
  double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
  tim=ru.ru_stime;        
  t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;        
  return t; 
} 

void difftime(double start, double end, unsigned int* time) {
  double secs = end - start;
  time[0] = (unsigned int)secs/3600;
  secs -= time[0]*3600;
  time[1] = (unsigned int)secs/60;
  secs -= time[1]*60;
  time[2] = secs;
}

