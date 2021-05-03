/*--------------------------------------------------------------------------*/
/*---------------------------- File EDPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the EDPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Niccolò Iardella, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "EDPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void EDPSolver::initialize( int k, ThermalUnitDPSolver * s ) {
 solver = s;
 auto & time_horizon = solver->time_horizon;

 h = k;
 kMax = time_horizon;

 int coeffsize = time_horizon * time_horizon + h * h - 2 * h * time_horizon;
 if( coeffsize != coeffs.size() ) {
  coeffs.resize( coeffsize );

  int msize = coeffsize + time_horizon - h;
  m.resize( msize );

  v.resize( time_horizon );
  pos.resize( time_horizon );
  unc_p.resize( time_horizon );
  con_p.resize( time_horizon );
 }
}

/*--------------------------------------------------------------------------*/

void EDPSolver::compute_costs( std::vector< double > & costs ) {

 // Scalar values
 auto & time_horizon = solver->time_horizon;
 auto & init_up_down_time = solver->init_up_down_time;
 auto & initial_power = solver->initial_power;

 // Power vectors
 auto & min_power = solver->min_power;
 auto & max_power = solver->max_power;
 auto & delta_ramp_up = solver->delta_ramp_up;
 auto & delta_ramp_down = solver->delta_ramp_down;
 auto & bound_on = solver->bound_on;
 auto & bound_down = solver->bound_down;

 // Coefficients of the objective function
 auto & quad_term = solver->quad_term;
 auto & linear_term = solver->linear_term;

 int k = h;

 coeffs[ 0 ].alfa = quad_term[ k ];
 coeffs[ 0 ].beta = linear_term[ k ];
 coeffs[ 0 ].gamma = 0; // TODO: Why 0?
 int coeffcnt = 1; // Next free position in coeffs[]
 v[ k ] = 0;       // Because for k = h the number of pieces is 1

 /*
  * Initialize the vector m containing the endpoints of the pieces.
  * At first, it contains the two endpoints of the individual piece.
  * At startup, power can't exceed the bound-on value /bar{l}_k.
  */

 if( h == 0 && init_up_down_time > 0 ) {
  m[ 0 ] = std::max( min_power[ k ], initial_power - delta_ramp_down[ k ] );
  m[ 1 ] = std::min( max_power[ k ], initial_power + delta_ramp_up[ k ] );
 } else {
  m[ 0 ] = min_power[ k ];
  m[ 1 ] = std::min( bound_on[ k ], max_power[ k ] ); // \bar{l}_k;
 }
 int mcnt = 2; // Next free position in m[]

 /*
  * Initialize the vector pos containing the initial indices of the pieces.
  */

 pos[ k ].begm = 0;
 pos[ k ].begt = 0;

 /*
  * Initialize the vector of unconstrained power values.
  * Unconstrained means that power values are not constrained by bound_down[k].
  */
 {
  // tmp is p^{*}_{hk}
  double tmp = -coeffs[ 0 ].beta / ( 2 * coeffs[ 0 ].alfa );
  if( tmp < m[ 0 ] ) {
   unc_p[ k ] = m[ 0 ];
  } else if( tmp > m[ 1 ] ) {
   unc_p[ k ] = m[ 1 ];
  } else {
   unc_p[ k ] = tmp;
  }
 }

 /*
  * Initialize the vector of constrained power values, that will be
  * computed at each iteration.
  * Constrained means that they must be <= bound_down[k].
  */
 if( k < time_horizon - 1 && unc_p[ k ] > bound_down[ k + 1 ] ) {
  con_p[ k ] = bound_down[ k + 1 ];
 } else {
  con_p[ k ] = unc_p[ k ];
 }

 costs[ k ] = coeffs[ 0 ].alfa * con_p[ k ] * con_p[ k ] +
              coeffs[ 0 ].beta * con_p[ k ];

 // Outermost loop
 for( k = h + 1; k < kMax; ++k ) {

  /*
   * Building pieces: \bar{m}_0 is the first endpoint of the first piece of
   * the z_{hk}(\bar{p}) objective function. Such endpoint will be saved in
   * the m vector.
   */

  pos[ k ].begm = mcnt;
  pos[ k ].begt = coeffcnt;


  if( min_power[ k ] > m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ] ) {
   m[ mcnt ] = min_power[ k ];
  } else {
   m[ mcnt ] = m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ];
  }

  double p_bar = m[ mcnt ]; // \bar{m}_0
  int v_bar = 0;            // After the case 3 will contain v[k]

  /*
   * Compute q, the index of the piece where p^*(\bar{p}) belongs.
   */

  double pstar; // p^*(\bar{p})

  if( p_bar < unc_p[ k - 1 ] ) {
   pstar = p_bar + delta_ramp_down[ k - 1 ];
   if( pstar > unc_p[ k - 1 ] ) {
    pstar = unc_p[ k - 1 ];
   }
  } else {
   pstar = p_bar - delta_ramp_up[ k - 1 ];
   if( pstar < unc_p[ k - 1 ] ) {
    pstar = unc_p[ k - 1 ];
   }
  }

  int qm = pos[ k - 1 ].begm;
  while( pstar >= m[ qm + 1 ] && qm < pos[ k ].begm - 2 ) {
   ++qm;
  }

  int q = qm - pos[ k - 1 ].begm + pos[ k - 1 ].begt;

  /*
   * Compute the last endpoint of the piece, \bar{u}.
   */

  double u_bar = std::min( max_power[ k ],
                           m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] );
  // if( max_power[ k ] < m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] ) {
  //  u_bar = max_power[ k ];
  // } else {
  //  u_bar = m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ];
  // }
  ++mcnt;


  int firstTime = 1;

  // CASE 1
  while( unc_p[ k - 1 ] > p_bar + delta_ramp_down[ k - 1 ] + eps ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta =
    linear_term[ k ] +
    coeffs[ q ].beta +
    2 * delta_ramp_down[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_down[ k - 1 ] * delta_ramp_down[ k - 1 ] +
    coeffs[ q ].beta * delta_ramp_down[ k - 1 ];

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible.
    */

   if( m[ qm + 1 ] - delta_ramp_down[ k - 1 ] <
       unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ] - eps ) {
    p_bar = m[ qm + 1 ] - delta_ramp_down[ k - 1 ];
    ++q;
    ++qm;
   } else {
    p_bar = unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ];
   }
   if( p_bar > u_bar ) {
    p_bar = u_bar;
   }
   ++v_bar;
   m[ mcnt++ ] = p_bar;

   /*
    * Compute unc_p, unconstrained optimal value for z_{hk}.
    */

   if( firstTime &&
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
    if( unc_p[ k ] < m[ mcnt - 2 ] ) {
     unc_p[ k ] = m[ mcnt - 2 ];
    }
    firstTime = 0;
   }

   ++coeffcnt;
  }

  // CASE 2
  if( unc_p[ k - 1 ] >= p_bar - delta_ramp_up[ k - 1 ] ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ];
   coeffs[ coeffcnt ].beta = linear_term[ k ];
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].alfa * unc_p[ k - 1 ] * unc_p[ k - 1 ] +
    coeffs[ q ].beta * unc_p[ k - 1 ] +
    coeffs[ q ].gamma;

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible.
    */

   if( ( unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ] ) < u_bar ) {
    p_bar = unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ];
   } else {
    p_bar = u_bar;
   }
   ++v_bar;
   m[ mcnt++ ] = p_bar;

   if( firstTime &&
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
    if( unc_p[ k ] < m[ mcnt - 2 ] ) {
     unc_p[ k ] = m[ mcnt - 2 ];
    }
    firstTime = 0;
   }

   ++coeffcnt;
  }


  // CASE 3
  while( p_bar < u_bar ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta =
    linear_term[ k ] + coeffs[ q ].beta -
    2 * delta_ramp_up[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_up[ k - 1 ] * delta_ramp_up[ k - 1 ] -
    coeffs[ q ].beta * delta_ramp_up[ k - 1 ];

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - \bar{p} stays admissible.
    */

   if( m[ qm + 1 ] + delta_ramp_up[ k - 1 ] < u_bar ) {
    p_bar = m[ qm + 1 ] + delta_ramp_up[ k - 1 ];
   } else {
    p_bar = u_bar;
   }
   ++v_bar;
   m[ mcnt++ ] = p_bar;
   ++q;
   ++qm;

   if( firstTime &&
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
    if( unc_p[ k ] < m[ mcnt - 2 ] ) {
     unc_p[ k ] = m[ mcnt - 2 ];
    }
    firstTime = 0;
   }

   ++coeffcnt;
  }
  // End of the tree cases

  v[ k ] = v_bar - 1;


  if( firstTime ) {
   // Function is strictly decreasing
   unc_p[ k ] = u_bar;
  }

  /*
   * Compute con_p[k], constrained optimal value for the entire function.
   */

  if( k < time_horizon - 1 && unc_p[ k ] > bound_down[ k + 1 ] ) {
   con_p[ k ] = bound_down[ k + 1 ];
  } else {
   con_p[ k ] = unc_p[ k ];
  }

  /*
   * Compute the cost for the node (h,k) in costs[].
   */

  qm = pos[ k ].begm;
  while( con_p[ k ] > m[ qm + 1 ] && m[ qm + 1 ] != 0 ) {
   ++qm;
  }
  q = qm - pos[ k ].begm + pos[ k ].begt;

  costs[ k ] =
   coeffs[ q ].alfa * con_p[ k ] * con_p[ k ] +
   coeffs[ q ].beta * con_p[ k ] +
   coeffs[ q ].gamma;
 }
}

/*--------------------------------------------------------------------------*/

void EDPSolver::compute_power_variables( int k, std::vector< double > & p ) {

 int t;
 auto & delta_ramp_up = solver->delta_ramp_up;
 auto & delta_ramp_down = solver->delta_ramp_down;

 /*
  * Compute optimal power values p_h[h], p_h[h+1], ..., p_h[k-1].
  */

 p[ k ] = con_p[ k ];
 for( t = k - 1; t >= h; --t ) {

  /*
   * Project unconstrained optimal value unc_p[t] on the interval:
   * [ p[t+1] - delta_ramp_up[t], p[t+1] + delta_ramp_down[t] ]
   *
   * If the unconstrained optimal value is on the left of the interval,
   * then the optimal power value is the left endpoint of the function.
   *
   * If the unconstrained optimal value is inside the interval,
   * then the optimal power value is exactly the unconstrained optimal value.
   *
   * If the unconstrained optimal value is on the right of the interval,
   * then the optimal power value is the right endpoint of the function.
   */

  if( unc_p[ t ] < p[ t + 1 ] - delta_ramp_up[ t ] ) {
   p[ t ] = p[ t + 1 ] - delta_ramp_up[ t ];
  } else if( unc_p[ t ] <= p[ t + 1 ] + delta_ramp_down[ t ] ) {
   p[ t ] = unc_p[ t ];
  } else {
   p[ t ] = p[ t + 1 ] + delta_ramp_down[ t ];
  }
 }
}

/*--------------------------------------------------------------------------*/
