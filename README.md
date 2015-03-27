# IMPRECISE HIDDEN MARKOV MODEL INFERENCE LIBRARY #

This C++ library implements inference algorithms for hidden Markov
models with imprecisely specified (i.e. set-valued) parameters. 


## OVERVIEW ##

This set of C++ files implements algorithms for (direct) inference in imprecise Hidden Markov Models (iHMMs). The currently implemented algorithms are:

- upper and lower likelihood probability computation (likelihood)
- upper and lower predictive probability (filtering)
- maximin and maximax Viterbi-like most likely state sequence estimation (viterbi)

## INSTALLATION ##

There is no installation. Simply run 'make' to compile the code from a terminal. The executable files are stored in the ./bin/ directory (assuming you are in the library's directory).

### FILES

The following executables are stored in the directory bin:

- likelihood,  likelihood computation [homgeneous models]
- viterbi, maximin and maximax state sequence estimation [homogeneous models]
- filter, filtering (calculate next hidden state probability) [homogeneous models]
- filter2, filtering (calculate next hidden state probability) [non-homogeneous models]

## DATA FORMAT ##

The algorithms accept iHHMs specified thorugh interval-valued probabilities. Most algorithms have only variants for homogeneous. The intervals are specified through tables of lower and upper transition, emission and prior probabilities in the following syntax (homogeneous models):

    # MODEL INPUT FILE TEMPLATE
    # lines starting with '#' are ignored
    # number of states
    N
    # number of outcomes
    M
    # lower transition probabilities
    p1 p2 ... pN pN+1 ... pN^2
    # upper transition probabilities
    p1 p2 ... pN pN+1 ... pN^2
    # lower emission probabilities
    p1 p2 ... pN pN+1 ... pN^2
    # upper emission probabilities
    p1 p2 ... pM pM+1 ... pN*M
    # lower prior probabilities
    p1 p2 ... pN
    # upper prior probabilities
    p1 p2 ... pN
    # blank or comment line at the end

In the above file, N denotes the number of values of hidden variables and M the number of values of manifest bariables. For transition probabilities, p[N*j+i] denotes Pr(Q'=i|Q=j), so e.g. pN+1 = Pr(Q'=1|Q=2). For emission probabilities p[N*i+j] = Pr(O=j|Q=i), and for prior probabilities p[i] = Pr(Q=i). The syntax for non-homogenous models is similar, with the first line denoting the number of steps, and each row describing distributions containing the probabilities for all time steps separated by spaces.

See files `weather.model` and `cweather.model` for examples of homogeneous (precise) HMMs and iHMMs. Non-homogeneous models are specified in a syntax similar to the above, except that the first line contains the horizon (number of steps), the number of states and symbols might change in each step (hence need to be specified), as well as the transition and emission probability itnervals. See `weather2.model` for an non-homogenous version of `weather.model` (with five time steps).

The data or observations (i.e., sequence of values for manifest variables) must be provided in a separate file using the syntax:

    # DATA FILE TEMPLATE 
    1st observation sequence
    2nd observation sequence
    3rd observation sequence
    ...
    Dth observation sequence
    # blank or comment line

See the file `weather.data` for an example.


## USAGE ##

Running an executable with no parameters outputs its usage.

## AUTHOR ##

(C) 2015 Denis D. Maua. You can contact me through denis.maua@gmail.com.

## LICENSE ##

This code is distributed under GPL 2 license. See file LICENSE.
