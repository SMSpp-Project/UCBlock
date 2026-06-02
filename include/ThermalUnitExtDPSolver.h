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
 * \author Claude Opus 4.7 \n
 *         Anthropic
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

 /// stages of the computation performed by compute()
 /** Each value records how far along the forward pipeline we have already
  * progressed for the current values of the ThermalUnitBlock parameters.
  * Once the Modification queue is drained, compute() only re-runs the
  * stages that have been invalidated (i.e., only those strictly after
  * the current stage value). Modifications that change input data reset
  * stage to start, forcing a full re-execution. */
 enum stage_value {
  start = 0 ,    ///< nothing computed (or parameters have changed)
  loaded_OK = 1 ,///< parameters loaded from ThermalUnitBlock (reserved)
  dp_OK = 2 ,    ///< forward DP has been run; f_best_cost is valid
  sol_OK = 3     ///< P[] and U[] have been assembled by backtracking
 };

/*--------------------------------------------------------------------------*/
 /// one quadratic piece alfa p^2 + beta p + gamma on the interval
 /// [left, right] of the power variable p
 /** This is the atomic building block of the piecewise quadratic value
  * functions used by the DP. The three coefficients encode the full
  * univariate quadratic alfa p^2 + beta p + gamma; the two endpoints
  * [left, right] delimit the portion of the power axis on which this
  * expression is meaningful. Neighbouring pieces in a PQFun share their
  * adjacent endpoints (right of piece i == left of piece i+1), so each
  * endpoint is stored redundantly but the representation is local and
  * self-contained, which makes insertion/splitting cheap during
  * sliding_min(). */

 struct PieceQuad {
  double alfa;   ///< coefficient of p^2 (>= 0 for a convex piece)
  double beta;   ///< coefficient of p
  double gamma;  ///< additive constant (absorbs accumulated path cost)
  double left;   ///< left endpoint of the piece (inclusive)
  double right;  ///< right endpoint of the piece (inclusive)
 };

 /// a convex piecewise quadratic function on the power axis
 /** Stored as a vector of PieceQuad sorted by left endpoint and covering a
  * contiguous sub-interval of the power axis with no gaps and no overlaps
  * (except at shared endpoints). The DP invariants are:
  *
  *  1. pieces[ i ].right == pieces[ i+1 ].left for all i (continuous support);
  *  2. alfa >= 0 on every piece (each piece is itself convex);
  *  3. the values at shared endpoints agree (function is continuous);
  *  4. the subgradient is monotonically non-decreasing across endpoints,
  *     so the *overall* piecewise function is convex.
  *
  * The empty vector represents the +INFinity function, used to signal an
  * infeasible or not-yet-reached DP state (e.g., a tau value that has no
  * predecessor chain, or an interval that is empty after clamping).
  *
  * All the "shape" operations implemented on PQFun (sliding_min, add_quadratic,
  * clamp_domain, ...) preserve invariants (1)-(4) by construction. */

 using PQFun = std::vector< PieceQuad >;

/*--------------------------------------------------------------------------*/
 /// summary of a single F^tau_t function, kept alongside the PQFun itself
 /** Produced by run_DP() right after each F^tau_t has been built. It is
  * used both to finalise the best cost at the end of the horizon and to
  * support backtracking without re-scanning the piecewise representation.
  *
  *  - min_val : the minimum value of F^tau_t(p) over the whole domain
  *              [min_power[t], max_power[t]] (or [min_power[t], bound_on[t]]
  *              when tau == 1, since a just-restarted unit is further
  *              constrained by the start-up ramp limit). TUEDPINF if the
  *              slot is infeasible (empty F).
  *  - argmin_p : a minimiser p_star of F^tau_t on the same domain. This is
  *              the "p*_{t}" referenced by eq. (16) of Wuijts et al. (2021):
  *              when walking the optimal path back from time t to time t-1,
  *              the power at t-1 is clamped to [p - Delta+, p + Delta-]
  *              around the current power p, and the "middle" branch of
  *              the three-case formula picks exactly p_star_prev.
  *
  * Predecessors are implicit and dictated by the forward DP recurrence:
  * for tau > 1 the predecessor of (t, tau) is (t-1, tau-1); for tau == 1
  * the predecessor is the off-side state at t-1 (specifically the
  * "ready-to-restart" state encoded by c_off_ready[t-1]). */

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
 /// take a spare PQFun from the pool (empty, capacity retained) or a new one
 PQFun pool_take( void ) {
  if( m_pqpool.empty() )
   return( PQFun{} );
  PQFun f = std::move( m_pqpool.back() );
  m_pqpool.pop_back();
  f.clear();  // keeps capacity
  return( f );
  }

 /// return a PQFun's storage to the pool for later reuse
 void pool_give( PQFun & f ) {
  f.clear();  // keeps capacity
  m_pqpool.push_back( std::move( f ) );
  }

 /** The result is written into \p out (cleared first), reusing its capacity;
  * an internal scratch buffer (m_raw) is reused across calls too, so the hot
  * path performs no per-call allocation. Non-static for that reason. */
 void sliding_min( const PQFun & F ,
                   double ramp_up , double ramp_down ,
                   double lo , double hi , PQFun & out );

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

 char stage;                  ///< computation stage (see stage_value)
 bool f_solved;               ///< true iff f_best_cost encodes a feasible
                              ///< schedule (finite cost); set by run_DP()
 double f_best_cost;          ///< optimal objective value after run_DP()

 // -- ON-side value functions ------------------------------------------ //

 /// surviving F^tau_t functions, in sparse parallel-vector layout
 /** f_F[ t ] and f_tau[ t ] have the same length and are sorted by tau
  * ascending: f_F[ t ][ i ] is the piecewise convex quadratic value
  * function for the (t, f_tau[t][i]) DP state. "Surviving" means two
  * things at once:
  *
  *  - *reachable*: only tau values that can be obtained by some legal
  *    sequence of DP transitions from the initial state are present;
  *    unreachable tau values simply have no entry (as opposed to an
  *    "+INF placeholder"). This collapses the storage from O(n * K_max)
  *    down to O(n * |kept|), which is what the paper relies on for
  *    practical efficiency;
  *
  *  - *relevant*: among the tau >= M_up entries, those pointwise dominated
  *    by another tau' >= M_up entry have been pruned by the pairwise
  *    RRF+ check (Wuijts et al. 2021, Prop. 6.1 / is_dominated_by); a
  *    dominated function cannot appear in any optimal schedule since its
  *    downstream propagation through sliding_min is dominated too.
  *
  * f_tau[t] is kept strictly increasing, so membership can be tested with
  * std::lower_bound in O(log |kept|) during backtracking. */
 std::vector< std::vector< PQFun > > f_F;
 std::vector< std::vector< Index > > f_tau;

 /// per-slot summary (min_val, argmin_p) parallel to f_F[t] / f_tau[t]
 /** Built at the same time as f_F[t] so we do not have to re-scan the
  * PQFun during backtracking. f_on[t][i] summarises f_F[t][i] over the
  * full domain of that function. */
 std::vector< std::vector< OnSlot > > f_on;

 // -- allocation pooling for run_DP() ---------------------------------- //
 // run_DP() rebuilds the sparse ON-side state from scratch at every
 // re-solve. To avoid the malloc/free churn of the many short-lived PQFun
 // (one per surviving entry per time step, from sliding_min) and of the
 // per-step list buffers, we recycle storage across calls:
 //  - m_pqpool holds spare PQFun buffers (capacity retained); the per-step
 //    functions are taken from it and the previous solve's f_F[t] are
 //    drained back into it at the next reset;
 //  - m_new_F/m_new_tau/m_new_on are the per-step build buffers, swapped
 //    into f_F[t]/f_tau[t]/f_on[t] instead of freshly allocated;
 //  - m_raw is sliding_min()'s internal scratch.
 std::vector< PQFun >  m_pqpool;
 std::vector< PQFun >  m_new_F;
 std::vector< Index >  m_new_tau;
 std::vector< OnSlot > m_new_on;
 PQFun                 m_raw;

 // -- OFF-side scalars ------------------------------------------------- //

 /// c_off_ready[ t ] : min cost of a schedule that is off at time t AND
 /// has been off for at least min_down_time consecutive instants (so it
 /// is legal to restart at time t + 1). This is the state from which
 /// F^1_{t+1} is built when c_off_ready[t] is finite. Updated at each t
 /// as the min of three sources: (1) stay from c_off_ready[t-1];
 /// (2) fresh shutdown happened at end of time t - min_down_time;
 /// (3) initial-off trail, when init_up_down_time <= 0 and the number of
 /// consecutive off instants since before the horizon (i.e.,
 /// |init_up_down_time| + t + 1) has just reached min_down_time.
 std::vector< double > c_off_ready;

 /// c_off_any[ t ] : min cost of a schedule that is off at time t,
 /// regardless of how long it has been off. Used only at the end of the
 /// horizon, where a trailing off period is not charged any further cost
 /// and the min-down-time of the *current* off period is not of our
 /// concern (it will be paid in the next planning horizon, if any).
 /// Updated as min( c_off_any[t-1], v_shutdown[t-1] ).
 std::vector< double > c_off_any;

 /// v_shutdown[ h ] : cost of reaching the "long shutdown arc" at the
 /// end of time h, i.e., min_{ tau >= M_up , p in [P, SD_{h+1}] } of
 /// F^tau_h(p), plus the constant shutdown cost (currently 0 since
 /// ThermalUnitBlock does not expose a shutdown cost). Only well defined
 /// for h < time_horizon - 1; the entry for the last instant is
 /// unused (the unit is free to stay on until the horizon end without
 /// entering a shutdown trajectory). Arriving at v_shutdown[h] makes the
 /// unit "ready to restart" at time h + M_down via the long arc. */
 std::vector< double > v_shutdown;

 /// optimal (tau, p) that achieves v_shutdown[h]: needed to backtrack
 /// through the long shutdown arc. If v_shutdown[h] is +INF, the two
 /// fields are undefined.
 std::vector< Index  > v_shutdown_tau;
 std::vector< double > v_shutdown_p;

 /// "origin time" of c_off_ready[ t ]: records the time h whose shutdown
 /// produced c_off_ready[t] via the long shutdown arc, so that
 /// backtracking knows where to jump back to. Value is -1 if:
 /// (a) c_off_ready[t] is +INF (never reached), or (b) the "ready" state
 /// was inherited from the initial off trail (init_up_down_time <= 0 and
 /// enough off time had accumulated before the horizon), in which case
 /// no in-horizon shutdown is involved and backtracking stops at t = 0.
 std::vector< int > f_ready_pred;

 /// analogous to f_ready_pred[ t ] but for c_off_any[ t ]: time of the
 /// last shutdown that contributed, or -1 when no in-horizon shutdown
 /// has occurred (init_up_down_time <= 0 and the unit has been off from
 /// the start of the horizon).
 std::vector< int > f_any_pred;

 // ---- output ---------------------------------------------------------- //

 /// power profile of the optimal schedule, as filled by build_solution()
 std::vector< double > P;
 /// commitment profile of the optimal schedule, as filled by build_solution()
 std::vector< bool >   U;

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitExtDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitExtDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File ThermalUnitExtDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
