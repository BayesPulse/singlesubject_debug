/// 
/// 
/// birthdeath.c 
/// 
/// ----------------------------------------------------------------------------
/// 
/// Global variable definitions:
///      fitstart - The first time in hours that a pulse may occur
///      fitend   - The last time in hours that a pulse may occur
///      mmm      - Order statistic used for distribution of pulse locations.
///                 This is inputted by the user and is typically 3
/// 
/// ----------------------------------------------------------------------------
/// 
/// Subroutines: 
///      birth_death
///      mean_contribution
///      *mean_concentration
///      likelihood
///      *calc_death_rate
/// 

#include "deconvolution_main.h"
#include <float.h>

// Global variable definitions 
extern double fitstart;
extern double fitend;
extern int mmm;

/// 
/// birth_death: 
///      this runs the birth-death step of the BDMCMC process
/// 
///      ARGUMENTS: 
///          Node_type *list     - this is the current list of pulses that exist
///          double **ts         - this is the matrix of observed data (a column
///                                of times and a column of log(concentration) 
///          Common_parms *parms - the current values of the common parameters
///          int N               - the number of observations in **ts
///          double *likeli      - the current value of the likelihood
///          unsigned long *seed - seed values needed for the randon number
///                                generator
///          int iter            - which iteration are we on
/// 
///      RETURNS:   
///          None                - all updates are made internally
/// 
/// /*{{{*/
///
///
///-----------------------------------------------------------------------------
///
/// Variable definitions
///   i,j,k              : generic counters
///   remove             : if a death occurs, remove is a draw from the
///                        multinomial distribution and represents which pulse
///                        we will remove
///   num_node           : number of pulses
///   flag               : if a birth occurs, we must draw parameters for this
///                        new pulse; flag will be set to 0 until a positive
///                        mass and a positive width are drawn
///   aaa                : this counter ensures we do not run the birthdeath
///                        loop too many times
///   max_num_node=60    : this variable is set to 60 and ensures a death if we
///                        exceed 60 pulses
///   S                  : the sum of exponential draws; the birthdeath loop
///                        stops when S exceeds T
///   Birth_rate         : a constant birth rate needed for the process
///   T=1                : stop birthdeath loop when S exceeds T=1
///   full_likelihood    : current likelihood before a birth or death occurs
///   max                : part of the calculation for Death_rate
///   Death_rate         : total death rate
///   *death_rate        : vector of death rates for each pulse
///   position           : if a birth occurs, the pulse's location is position,
///                        drawn from a uniform distribution
///   *partial_likelihood: a vector of likelihoods where the ith element is the
///                        likelihood if pulse i were removed
///   *tmp               : used when drawing new pulse's mass and width
///   *node              : counter through pulses
///   *new_node          : if a birth occurs, we give the new pulse this name
///
// -----------------------------------------------------------------------------
//
// Subroutines used
//      rnorm             - found in randgen.c; draws from the normal
//                          distribution
//      rmultinomial      - found in randgen.c; draws from the multinomial
//                          distribution
//      kiss              - found in randgen.c; draws from U(0,1) distribution
//      rexp              - found in randgen.c; draws from the exponential
//                          distribution
//      runif_atob        - found in randgen.c; draws from U(a,b) distribution
//      mean_contribution - found in this file; evaluates mean_contrib for a
//                          pulse with specific parameters
//      likelihood        - found in this file; calculates and returns
//                          likelihood under current parameters 
//                          Node_type *initialize_node(void)
//      *initialize_node  - found in hash.c; allocates memory for a new pulse
//                          and creates it;
//      insert_node       - found in hash.c; inserts a newly created node in the
//                          appropriate spot in the linked list
//      delete_node       - found in hash.c; eliminates a node and links the
//                          neighbors
//      *calc_death_rate  - found in this file; creates a vector where the ith
//                          element represents the death rate for pulse i
//
///
void birth_death(Node_type *list, double **ts, Common_parms *parms, int N, 
                 double *likeli, unsigned long *seed, int iter) 
{

  // Declarations
  int i;
  int j; 
  int k; 
  int remove; 
  int num_node; 
  int flag; 
  int aaa = 0;  // BD iteration counter
  int max_num_node = 60;
  double S = 0.0; // Current progress in virtual time
  double Birth_rate;
  double T = 1.;
  double full_likelihood;
  double max;
  double Death_rate;
  double *death_rate;
  double position;
  double *partial_likelihood;
  double *tmp;
  Node_type *node;
  Node_type *new_node;

  double rnorm (double,double,unsigned long *); // normal RNG
  void print_list(Node_type *); // prints pulses linklist
  long rmultinomial (double *,long,unsigned long *); // multinomial RNG
  double kiss (unsigned long *); // uniform RNG
  double rexp (double, unsigned long *); // exponential RNG
  double runif_atob (unsigned long *,double,double); // unif(a,b) RNG
  // Calculates concentration at each timepoint for each pulse individually
  void mean_contribution (Node_type *,double **,Common_parms *,int); 
  // Calculates likelihood for given parms, option: exclude one with node_out
  double likelihood (Node_type *list,double **ts,Common_parms *parms,int N,
                     Node_type *node_out);
  // Create new pulse and associated parameters
  Node_type *initialize_node (void);
  // Insert new pulse into linklist in appropriate order
  void insert_node (Node_type *new_node,Node_type *list);
  // Delete pulse
  void delete_node (Node_type *node,Node_type *list);
  // Calculate death rate for each pulse
  double *calc_death_rate (Node_type *, int, double *, double, double, double);


  // Set birth rate to poisson prior lambda
  Birth_rate = parms->nprior;

  tmp = (double *)calloc(2, sizeof(double));

  // Save Likelihood
  full_likelihood = *likeli;

  // Go until the loop is broken
  while (1) {

    // This counter keeps this loop from running too many times
    aaa++;

    // Count number of pulses
    num_node = 0;
    node = list->succ;
    while (node != NULL) {
      num_node++;
      node = node->succ;
    }

    // Allocate memory for partial likelihood vector
    partial_likelihood = (double *)calloc(num_node,sizeof(double));

    // Calculate the likelihood if pulse i is removed
    i = 0;
    node = list->succ;
    while (node != NULL) {
      partial_likelihood[i] = likelihood(list, ts, parms, N, node);
      i++;
      node = node->succ;
    }


    // Calculate death rate for each component 
    death_rate = NULL;
    death_rate = calc_death_rate(list, num_node, partial_likelihood,
                                 full_likelihood, Birth_rate, parms->nprior);

    // Compute total death rate D
    //   The next portion computes D = sum(d_i); This is a little more
    //   complicated than summing them up because of precision issues.
    if (death_rate != NULL) {
      Death_rate = death_rate[0];

      for (i = 1; i < num_node; i++) {
        max = (Death_rate > death_rate[i]) ? Death_rate : death_rate[i];
        Death_rate = log(exp(Death_rate - max) + 
                         exp(death_rate[i] - max)) + max;
      }

      for (i = 0; i < num_node; i++) {
        death_rate[i] -= Death_rate;
        death_rate[i]  = exp(death_rate[i]);
      }

      for (i = 1; i < num_node; i++) {
        death_rate[i] += death_rate[i-1];
      }

      if (Death_rate > 500 ) {
        Death_rate = 1e300;
      } else {
        Death_rate = exp(Death_rate); 
      }

    } else { Death_rate = 0; }

    free(partial_likelihood);

    if (num_node <= 1) Death_rate = 0;

    // Jump in virtual time 
    S += rexp(Birth_rate + Death_rate, seed);

    // If S exceeds T or if we've run this too many times, move on
    if (S > T)  { 
      //printf("BD ran for %d iterations.\n", aaa);
      break;
    }
    if (aaa > 5000) {
      //printf("BD ran for %d iterations.\n", aaa);
      break;
    }

    // Calculate prob(birth)
    if (num_node <= 1) { // birth if 0 or 1 pulses
      max = 1.1;
    } else if (num_node >= max_num_node) { // death if too many exist
      max = -0.1;
    } else { // Otherwise, set max = B/(B+D) 
      max = Birth_rate / (Birth_rate + Death_rate); 
    }


    // If U < B/(B+D), a birth occurs 
    if (kiss(seed) < max) {  

      // Simulate a pulse position uniformly 
      position = runif_atob(seed, fitstart, fitend);

      node = list->succ;

      // Simulate random effects from the prior 

      for (j = 0; j < 2; j++){
        // lognormal prior on rand effects
        tmp[j] = exp(rnorm(parms->theta[j],
                           parms->re_sd[j], seed));
      } 

      // initialize a new node 
      new_node               = initialize_node();
      new_node->time         = position;
      new_node->theta[0]     = tmp[0];
      new_node->theta[1]     = tmp[1];
      new_node->mean_contrib = (double *)calloc(N, sizeof(double));
      mean_contribution(new_node, ts, parms, N);
      // insert the new node in correct place 
      insert_node(new_node, list);
      full_likelihood = likelihood(list, ts, parms, N, list);

    } else { // Otherwise, a death occurs 

      // Draw node to kill from multinomial
      remove = (int)rmultinomial(death_rate,(long)num_node,seed) + 1;
      node = list;
      for (i = 0; i < remove; i++) { node = node->succ; }
      delete_node(node,list);
      full_likelihood = likelihood(list, ts, parms, N, list);

    }

    free(death_rate); //NOTE: necessary?

  } // End of While Loop

  // Update likelihood before leaving Birthdeath.c
  *likeli = full_likelihood;

  free(death_rate);
  free(tmp);

}
/*}}}*/




/*! 
 * mean_contribution: 
 *      this updates a pulse's mean_contrib vector based on current values of
 *      parameters
 *
 *      Two version: Phi version and built-in erf version
 *
 *      ARGUMENTS: 
 *          Node_type *node     - what pulse are we updating
 *          double **ts         - this is the matrix of observed data (a column
 *                                of times and a column of log(concentration)
 *          Common_parms *parms - the current values of the common parameters
 *          int N               - the number of observations in **ts
 *
 *      RETURNS: 
 *          None                - all updates are made internally
 *
 */ /*{{{*/
/*!
 *
 * -----------------------------------------------------------------------------
 *
 * Variable definitions
 *      i         - generic counter
 *      x,y,z,tmp - these are all variables that are part of the arithmetic
 *                  needed in calculating the mean contribution
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines used
 *    --------erf         - this is a subroutine used to help in integration
 *
 */
// Phi() version
//void mean_contribution(Node_type *node, double **ts, 
//                       Common_parms *parms, int N){
//    /* 
//     * compute the contribution to the mean hormone concentration from a pulse
//     */
//    int i;
//    double x,y,z,twopi=0.6931471805599453;
//    double Phi(double);
//    
//    z = node->theta[1]*twopi/parms->md[1];
//    y = twopi*(0.5*z/parms->md[1] + node->time/parms->md[1]);
//    z += node->time;
//
//    for (i=0;i<N;i++) {
//        x = Phi((ts[i][0]-z)/sqrt(node->theta[1])) - Phi(-z/sqrt(node->theta[1]));
//        node->mean_contrib[i] = node->theta[0]*x*exp(y - ts[i][0]*twopi/parms->md[1]);
//
//        if (!(node->mean_contrib[i] > -DBL_MAX && node->mean_contrib[i] < DBL_MAX)){
//             node->mean_contrib[i] = 0;
//        }
//        
//    } 
//
//}

// erf() version -- via code9uniform 
void mean_contribution(Node_type *node, double **ts, Common_parms *parms, int N)
{
  /* This function updates a pulse's mean_contrib vector based on inputted parms*/
  int i;
  double x, y, z, tmp; //, a, temp;
  double erf(double);

  z = node->theta[1]*0.6931472/parms->md[1];
  y = 0.6931472*(0.5*z/parms->md[1] + node->time/parms->md[1]);
  z += node->time;
  tmp = sqrt(2.*node->theta[1]);

  for (i=0;i<N;i++) {
    x = (ts[i][0] - z)/tmp;
    x = 0.5*(1.+erf(x));

    if (x == 0) {
      node->mean_contrib[i] = 0; 
    } else {
      node->mean_contrib[i] = 
        node->theta[0] * x * exp(y - ts[i][0] * 0.6931472 / parms->md[1]);
    }

  }
}
/*}}}*/


/*! 
 * phi: 
 *      Returns the normal distribution function (with mean mu, sd s) evaluated
 *      at y, parameters
 *
 *      ARGUMENTS: 
 *          double y            - value to evaluate at
 *          double mu           - distribution mean
 *          double s            - distribution standard deviation (variance?)
 *
 *      RETURNS: 
 *          *value*             - 
 *
 */ /*{{{*/
/*!
 *
 * -----------------------------------------------------------------------------
 *
 * Variable definitions
 *      x,t,z     - these are all variables that are part of the arithmetic
 *                  needed in calculating the mean contribution
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines used
 *
 */
//double phi(double y, double mu, double s){
//  /* Returns normal distribution function evaluated at y */
//  /* Returns standard normal distribution function evaluated at y */
//  /* Has absolute error of order 10^(-7) */
//  /* Uses approximation to erfc from Numerical Recipes in C */
//  double t, z, ans, x;
//
//  x = (y-mu)/(1.4142135623730951*s);
//  z = fabs(x);
//  t = 1.0/(1.0+.5*z);
//  ans = t * exp(-z * z-1.26551223 + t * 
//                (1.00002368 + t * 
//                 (0.37409196 + t * 
//                  (0.09678418 +  t * 
//                   (-0.18628806 + t * 
//                    (0.27886807 + t * 
//                     (-1.13520398 + t * 
//                      (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))); 
//  if (x >= 0){
//    return 1.0-0.5*ans;
//  } else {
//    return 0.5*ans;
//  }
//
//}
/*}}}*/




/*!
 * Phi
 *      standard normal version of phi
 */
/*{{{*/
//double Phi(double y) {
//  /* Returns standard normal distribution function evaluated at y */
//  /* Has absolute error of order 10^(-7) */
//  /* Uses approximation to erfc from Numerical Recipes in C */
//  double t, z, ans, x, inv_sqrt2 = 0.7071067811865475;
//  x = y*inv_sqrt2;
//  z = fabs(x);
//  t = 1.0/(1.0+.5*z);
//  ans = t * exp(-z * z-1.26551223 + t * 
//                (1.00002368 + t * 
//                 (0.37409196 + t * 
//                  (0.09678418 +  t * 
//                   (-0.18628806 + t * 
//                    (0.27886807 + t * 
//                     (-1.13520398 + t * 
//                      (1.48851587 + t * 
//                       (-0.82215223 + t * 0.17087277)))))))));
//
//  if (!(x < 0)){
//    return 1.0-0.5*ans;
//  } else {
//    return 0.5*ans;
//  }
//
//}
/*}}}*/




/*!
 * mean_concentration: 
 *      this takes each pulse's mean_contrib vector and sums across them
 *
 *      ARGUMENTS: 
 *          Node_type *list     - this is the current list of pulses that exist
 *          Common_parms *parms - the current values of the common parameters
 *          int N               - the number of observations in **ts
 *          Node_type *node_out - if we want, we can ignore a pulse
 *          double **ts         - this is the matrix of observed data (a column
 *                                of times and a column of log(concentration)
 *
 *      RETURNS: 
 *          x                   - the vector of sums
 *
 */
/*{{{*/
/*!
 *
 * -----------------------------------------------------------------------------
 *
 * Variable definitions
 *      i     : generic counter
 *      *x    : the vector of sums
 *      *node : the counter of pulses
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines used
 *      rnorm : found in randgen.c; draws from the normal distribution
 *
 */

double *mean_concentration(Node_type *list, Common_parms *parms, int N,
                           Node_type *node_out, double **ts){

  /* 
   * This function sums the mean_contrib vector across pulses 
   * The result is a vector of sums 
   */
  int i;
  double *x;
  Node_type *node;
  double rnorm(double, double, unsigned long *);

  x = (double *)calloc(N, sizeof(double));

  /* add the contribution to the mean from each pulse */
  node = list->succ;
  while (node != NULL) {
    if (node != node_out) {
      for (i=0; i<N; i++) {
        x[i] += node->mean_contrib[i];
      }
    }
    node = node->succ;
  }

  /* add the baseline contribution */
  for (i=0; i<N; i++) {
    x[i] += parms->md[0];
    x[i] = log(x[i]);
  }

  return x;

}
/*}}}*/




/*!
 * likelihood: 
 *      computes the current likelihood using the observed log-concentrations
 *      and mean concentration
 *
 *      ARGUMENTS: 
 *          Node_type *list     - this is the current list of pulses that exist
 *          double **ts         - this is the matrix of observed data (a column
 *                                of times and a column of log(concentration)
 *          Common_parms *parms - the current values of the common parameters
 *          int N               - the number of observations in **ts
 *          Node_type *node_out - if we want, we can ignore a pulse
 *
 *    RETURNS: 
 *          x                   - the likelihood computed based on the inputs
 *
 */ 
/*{{{*/
/*
 * 
 * -----------------------------------------------------------------------------
 *
 * VARIABLE DEFINITIONS
 *      i     - generic counter
 *      x=0   - likelihood to be calculated; initialized at zero
 *      *mean - vector of sums of mean_contribs calculated using
 *              mean_concentration
 *
 * -----------------------------------------------------------------------------
 *
 * SUBROUTINES USED
 *      mean_concentration - found in this file; computes a vector of sums of
 *                           mean_contribs
 *
 */

double likelihood(Node_type *list, double **ts, Common_parms *parms, int N,
                  Node_type *node_out){

  /* This function computes the likelihood under inputted parameters 
   * The output is a scalar */
  int i;
  double x=0,*mean;

  double *mean_concentration(Node_type *, Common_parms *, int, Node_type *,
                             double **);

  /* Sum across mean_contribs */
  mean = mean_concentration(list, parms, N, node_out, ts);

  for (i=0;i<N;i++) {
    x += (ts[i][1]-mean[i])*(ts[i][1]-mean[i]);
  }

  x /= (-2.0*parms->sigma);
  x += -0.5*N*(1.8378771+parms->lsigma);
  free(mean);

  return x;

}
/*}}}*/




/*
 * calc_death_rate: 
 *      calculates a vector of death rates, one for each existing pulse
 *    
 *      ARGUMENTS: 
 *          Node_type *list            - this is the current list of pulses that
 *                                       exist
 *          int num_node               - current number of pulses
 *          double *partial_likelihood - vector of likelihoods, where the ith
 *                                       element represents the likelihood with
 *                                       the ith pulse removed
 *          double full_likelihood     - value of the full likelihood
 *          double Birth_rate          - value of the birth rate
 *
 *      RETURNS:
 *          death_rate                 - a vector where the ith element
 *                                       represents the death rate of the ith
 *                                       pulse
 *
 */ /*{{{*/
/*
 *  
 * -----------------------------------------------------------------------------
 *
 * Variable definitions
 *      i,j                  - generic counters
 *      x                    - variable for death rate
 *      *death_rate          - vector of death rates
 *      coef_denom, coef_num - part of the calculation of death rate
 *      *node                - counter for pulses
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines used
 *      None
 *
 */

double *calc_death_rate(Node_type *list, int num_node, 
                        double *partial_likelihood, double full_likelihood, 
                        double Birth_rate, double r) {

  /* This function calculates the death rate vector */
  int i,j;
  double x,*death_rate;
  double coef_denom,coef_num;
  Node_type *node;

  if (num_node > 1) {
    death_rate = (double *)calloc(num_node, sizeof(double));
    node = list->succ;
    i = 0;

    /* Calculate the coefficient of the distribution of the taus conditional
       on number of pulses. In this portion, I have an extra num_node in the
       numerator */
    coef_num = 1;
    for (j=1;j<mmm;j++) coef_num *= j;
    coef_denom = mmm*num_node;
    for (j=1;j<mmm;j++) coef_denom *= (mmm*num_node+j);

    while (node != NULL) {

      if (mmm > 0) {
        /* This computes the portion of death rate due to the Poisson prior on
           N, the birth rate, and the likelihood */
        // num_node = number of pulses
        // r = prior on pulse count
        x = log(num_node*Birth_rate/r) + partial_likelihood[i] -
          full_likelihood;
      } else {
        x = log(Birth_rate/num_node) + partial_likelihood[i] -
          full_likelihood;
      }

      /*Now, compute the portion of the death rate due to distribution of the
        taus conditional on number of pulses. */
      if (mmm > 1) {
        if (i==0) {

          /*If we are on the first pulse*/
          x += (mmm-1)*log(fitend - fitstart) +
            (mmm-1)*log((node->succ->time - fitstart) /
                        ((node->succ->time - node->time) * 
                         (node->time - fitstart))) +
            log(coef_num/coef_denom);

        } else if (i>0) {

          /*If we are not on the first pulse */
          if (node->succ) {
            x += (mmm-1)*log(fitend - fitstart) +
              (mmm-1)*log((node->succ->time - node->pred->time) /
                          ((node->succ->time - node->time) * 
                           (node->time - node->pred->time))) +
              log(coef_num/coef_denom);

            /*If we are not on the first or last pulse */
          } else {
            x += (mmm-1)*log(fitend - fitstart) + 
              (mmm-1)*log((fitend - node->pred->time) / 
                          ((fitend - node->time) * 
                           (node->time - node->pred->time))) +
              log(coef_num/coef_denom);
          }

        }
      }

      // if pulse equal to NaN, then set to large value/guarantee death
      // necessary when starting value is < fitstart (orderstat part of
      // death calc causes this due to neg. value in log) 
      if (isnan(x)) {
        x = 1e300;
      }
      /*Save to death rate vector */
      death_rate[i] = x;

      /*Advance to next pulse*/
      node = node->succ;

      /*Advance counter 1*/
      i++;

    }

  } else {

    /*if we have 0 or 1 pulses*/
    if (num_node==0) {
      /*if we don't have any pulses, no death rate*/
      return NULL;
    } else {
      /*If we have 1 pulse, don't kill it*/
      death_rate = (double *)calloc(num_node,sizeof(double));
      death_rate[0] = -1e300;
    }
  }

  return death_rate;  

}
/*}}}*/




/*******************************************************************************
 * END OF FILE
 ******************************************************************************/

