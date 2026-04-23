/*--------------------------------------------------------------------------*/
/*---------------------- File ThermalUnitExtDPSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitExtDPSolver class, a Solver for the
 * ThermalUnitBlock (without primary and secondary reserve variables) that
 * solves the single-Unit Commitment (1UC) problem by a "hybrid" Dynamic
 * Programming scheme, broadly inspired by Wuijts, van den Akker and van den
 * Broek (Electric Power Systems Research, 2021) but with the off-state of
 * the commitment state-space collapsed to a single layer.
 *
 * In more detail:
 *
 * - The ON portion of the state-space is modelled as in the paper: at each
 *   time instant t and for each run-length tau (number of consecutive time
 *   instants for which the unit has been on, ending with t), a *convex
 *   piecewise quadratic function*
 *
 *                   F^tau_t : [P_min_t, P_max_t] -> R
 *
 *   stores the optimal cost of a schedule that is on at t with power p and
 *   whose current on-run has length exactly tau. These functions are built
 *   inductively by standard "ramp-constrained sliding minimum" (equivalently,
 *   the same one-step transition of Frangioni and Gentile, 2006) plus the
 *   addition of the cost f_t(p) = alfa_t p^2 + beta_t p + gamma_t. Transitions
 *   are arcs of length one between (on^{tau-1}, t-1) and (on^tau, t).
 *
 * - The OFF portion of the state-space collapses Wuijts's M_down counter
 *   into a *single* layer: at each time instant t, two scalars are kept,
 *
 *       c_off_ready(t) = optimal cost of a schedule that is off at t AND
 *                        has been off for at least M_down consecutive steps
 *                        (i.e., it is legal to restart at t+1);
 *
 *       c_off_any(t)   = optimal cost of a schedule that is off at t,
 *                        regardless of how long (used only at the end of
 *                        the horizon where an off trailing tail is free).
 *
 *   The min down-time constraint is enforced by a "long" shutdown arc that
 *   skips M_down instants ahead: a shutdown decided at the end of time h
 *   contributes to c_off_ready(h + M_down), not to c_off_ready(h + 1).
 *
 * - Transitions:
 *
 *     on  -> on  :  F^tau_t(p) = f_t(p) + min_{q in [p-Delta+, p+Delta-]}
 *                                F^{tau-1}_{t-1}(q)       (1 < tau <= K_max)
 *     off -> on  :  F^1_t(p)   = f_t(p) + SUC(t) + c_off_ready(t-1)
 *                                                           (p in [P, SU])
 *     on  -> off :  v_shutdown(h) = min_{p in [P, SD]} min_{tau >= M_up}
 *                                   F^tau_h(p) + c_stop
 *                   c_off_ready(h + M_down) <- min ..., v_shutdown(h)
 *     off -> off :  c_off_ready(t) = min( c_off_ready(t-1), v_shutdown(t - M_down) )
 *                   c_off_any(t)   = min( c_off_any(t-1),   v_shutdown(t - 1) )
 *
 * tau is in { 1, ..., K_max }, where K_max is bounded by max(t+1, init+t+1),
 * i.e., full RRF-style book-keeping without merging into an absorbing state;
 * all F^tau_t stay convex and the sliding minimum is computed in closed form
 * (cf. Wuijts eq. (16) or Frangioni-Gentile (2006)).
 *
 * This is the "multi-layer ON / single-layer OFF" hybrid discussed in
 * tandem with ThermalUnitDPSolver, which instead uses a single layer on
 * both sides and long arcs on both sides. Data loading, Modification
 * handling and Solution writing follow the same pattern as
 * ThermalUnitDPSolver (to which this class is otherwise unrelated).
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ThermalUnitExtDPSolver
 #define __ThermalUnitExtDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ThermalUnitExtDPSolver ----------------------*/
/*--------------------------------------------------------------------------*/
/// DP solver for 1UC with multi-layer ON / single-layer OFF graph

class ThermalUnitExtDPSolver : public Solver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 static constexpr auto TUEDPINF = Inf< double >();  ///< the INF value

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 ThermalUnitExtDPSolver( void ) {}

 ~ThermalUnitExtDPSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// sets the Block that the Solver has to solve
 void set_Block( Block * block ) override;

 /// solves the constructed problem
 int compute( bool changedvars = true ) override;

 /// tells whether a solution is available
 bool has_var_solution( void ) override { return( f_solved ); }

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc ) override;

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override { return( f_best_cost ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( f_best_cost ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override { return( f_best_cost ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// stage of the computation
 enum stage_value {
  start = 0 ,    ///< nothing computed
  loaded_OK = 1 ,///< parameters loaded
  dp_OK = 2 ,    ///< DP has been run
  sol_OK = 3     ///< solution P[], U[] has been assembled
 };

/*--------------------------------------------------------------------------*/
/// one quadratic piece alfa p^2 + beta p + gamma on [left, right]

 struct PieceQuad {
  double alfa;
  double beta;
  double gamma;
  double left;
  double right;
 };

 /// convex piecewise quadratic function, stored as sorted list of pieces
 /** Each element describes the piece alfa p^2 + beta p + gamma on
  * [left, right]. The pieces are sorted by left endpoint and cover a
  * contiguous interval with no gaps and no overlaps (except at shared
  * endpoints). Because the functions used in the DP are convex, alfa >= 0
  * on every piece and the piecewise junction is continuous. An empty
  * vector denotes the +INF function, i.e., an infeasible state. */

 using PQFun = std::vector< PieceQuad >;

/*--------------------------------------------------------------------------*/
/// backtracking record for a single F^tau_t "slot"
 /** After the forward DP each (t, tau) slot is summarised by:
  *
  *  - min_val : min over p of F^tau_t(p) on the full domain (INF if the slot
  *    is infeasible);
  *  - argmin_p : the corresponding unconstrained minimiser p_star of
  *    F^tau_t; this is the value used by the "middle" branch of the
  *    Wuijts sliding minimum while walking the optimal path back.
  *
  * Predecessors are implicit: for tau > 1 the predecessor is (t-1, tau-1);
  * for tau == 1 the predecessor is the off-side state at t-1. */

 struct OnSlot {
  double min_val;
  double argmin_p;
 };

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// read all the parameters from the ThermalUnitBlock
 void load_parameters( void );

 /// process the queue of Modifications
 void process_modifications( void );

 /// process one Modification; returns true if a full reload is required
 bool guts_of_process_modifications( const p_Mod mod );

 /// expand a per-instant vector as in ThermalUnitDPSolver::retrieve_term()
 void retrieve_term( std::vector< double > & out ,
                     const std::vector< double > & in ) const;

/*--------------------------------------------------------------------------*/

 /// run the full forward DP (ON and OFF layers together)
 void run_DP( void );

 /// reconstruct the commitment/power schedule by backtracking the DP
 void build_solution( void );

/*--------------------------------------------------------------------------*/
/*---------------- PIECEWISE-QUADRATIC FUNCTION HELPERS --------------------*/
/*--------------------------------------------------------------------------*/

 /// evaluate a quadratic piece coefficients (alfa, beta, gamma) at p
 static double eval_piece( const PieceQuad & pc , double p ) {
  return( pc.alfa * p * p + pc.beta * p + pc.gamma );
 }

 /// evaluate a piecewise quadratic function F at p; returns TUEDPINF if
 /// p is outside the domain
 static double eval( const PQFun & F , double p );

 /// argmin of a single quadratic piece on [left, right]
 /** Returns the value p in [pc.left, pc.right] that minimises
  * eval_piece(pc, .). */
 static double argmin_piece( const PieceQuad & pc );

 /// minimum of a PQFun over [lo, hi] intersected with its domain
 /** Returns { min_value, argmin_p }; if the intersection is empty, returns
  * { TUEDPINF, 0 }. */
 static std::pair< double , double > min_over(
  const PQFun & F , double lo , double hi );

 /// pointwise addition of the constant 'c' to F
 static void shift_by( PQFun & F , double c );

 /// pointwise addition of (alfa p^2 + beta p + gamma) to F
 static void add_quadratic( PQFun & F ,
                            double alfa , double beta , double gamma );

 /// restrict F to the domain [lo, hi] (in place)
 static void clamp_domain( PQFun & F , double lo , double hi );

 /// ramp-constrained sliding minimum of a convex piecewise quadratic F
 /** Computes G(p_t) = min_{q in [p_t - ramp_up, p_t + ramp_down]} F(q) on
  * the domain [lo, hi]. Assumes F convex (alfa >= 0 in every piece).
  * Implementation follows Wuijts et al. (2021) eq. (16)-(21), equivalent
  * to the three-case analysis of Frangioni and Gentile (2006). */
 static PQFun sliding_min( const PQFun & F ,
                           double ramp_up , double ramp_down ,
                           double lo , double hi );

 /// check whether F1 is pointwise >= F2 (up to tolerance eps) on all of
 /// dom( F1 ); if dom( F1 ) extends beyond dom( F2 ), returns false (at
 /// that point F2 is +INF which cannot dominate F1). Used by RRF+ to
 /// detect irrelevant F^tau functions.
 static bool is_dominated_by( const PQFun & F1 , const PQFun & F2 ,
                              double eps = 1e-9 );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // ---- data loaded from the ThermalUnitBlock --------------------------- //

 Index time_horizon;          ///< time horizon
 int init_up_down_time;       ///< initial up/down time (can be < 0)
 Index min_up_time;           ///< minimum up time
 Index min_down_time;         ///< minimum down time
 double initial_power;        ///< initial power
 Index t_init;                ///< first instant in which commitment is free

 std::vector< double > startup_costs;
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 std::vector< double > min_power;
 std::vector< double > max_power;
 std::vector< double > bound_on;
 std::vector< double > bound_down;

 std::vector< double > quad_term;
 std::vector< double > linear_term;
 std::vector< double > const_term;

 double eps{ 1e-10 };         ///< numerical tolerance

 // ---- DP state ------------------------------------------------------- //

 char stage;                  ///< what has been computed
 bool f_solved;               ///< set by run_DP()+build_solution()
 double f_best_cost;          ///< objective value of the current solution

 /// F^tau_t functions, sparse. f_F[ t ] and f_tau[ t ] are parallel
 /// vectors (same length) sorted by tau ascending: f_F[ t ][ i ] is the
 /// surviving function with tau == f_tau[ t ][ i ] at time t. Surviving
 /// means: reachable by some legal DP path AND not flagged as irrelevant
 /// by the pairwise-domination heuristic (RRF+).
 std::vector< std::vector< PQFun > > f_F;
 std::vector< std::vector< Index > > f_tau;

 /// per-slot summary for backtracking: parallel to f_F[ t ] / f_tau[ t ]
 std::vector< std::vector< OnSlot > > f_on;

 /// c_off_ready[ t ]: min cost of being off at t, "ready" (off >= M_down)
 std::vector< double > c_off_ready;

 /// c_off_any[ t ]: min cost of being off at t (any residual off duration)
 std::vector< double > c_off_any;

 /// v_shutdown[ t ] (in range [0, time_horizon - 1]): min cost of a
 /// schedule on through instant t with final power in [P, SD] plus the
 /// shutdown cost (cstop), i.e., the cost "arriving" at the long
 /// shutdown arc
 std::vector< double > v_shutdown;

 /// for each shutdown index, which (tau, p) achieved v_shutdown[ t ]
 /// (tau in [1, K_max[t]]; p in the shutdown range). Only used for
 /// backtracking. If v_shutdown[t] is +INF, both fields are undefined.
 std::vector< Index  > v_shutdown_tau;
 std::vector< double > v_shutdown_p;

 /// last shutdown time whose long arc produced c_off_ready[ t ], or -1 if
 /// c_off_ready[ t ] is +INF (or the "ready" came from the initial state)
 std::vector< int > f_ready_pred;

 /// last shutdown time contributing to c_off_any[ t ], or -1 if the unit
 /// has been off from before the horizon (no shutdown happened yet)
 std::vector< int > f_any_pred;

 // ---- output ---------------------------------------------------------- //

 std::vector< double > P;     ///< power values of the optimal schedule
 std::vector< bool >   U;     ///< commitment values of the optimal schedule

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitExtDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitExtDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File ThermalUnitExtDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
