/*--------------------------------------------------------------------------*/
/*---------------------- File ThermalUnitExtDPSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitExtDPSolver class, a Solver for the
 * ThermalUnitBlock that solves the single-Unit Commitment (1UC) problem by
 * a "hybrid" Dynamic Programming scheme, broadly inspired by Wuijts, van
 * den Akker and van den Broek (Electric Power Systems Research, 2021) but
 * with the off-state of the commitment state-space collapsed to a single
 * layer.
 *
 * In more detail:
 *
 * - The ON portion of the state-space is modelled as in the paper: at each
 *   time instant \f$ t \f$ and for each run-length \f$ \tau \f$ (number of
 *   consecutive time instants for which the unit has been on, ending with
 *   \f$ t \f$), a *convex piecewise quadratic function*
 *
 *   \f[ F^\tau_t : [ P^{min}_t , P^{max}_t ] \rightarrow \mathbb{R} \f]
 *
 *   stores the optimal cost of a schedule that is on at \f$ t \f$ with
 *   power \f$ p \f$ and whose current on-run has length exactly
 *   \f$ \tau \f$. These functions are built inductively by the standard
 *   "ramp-constrained sliding minimum" (equivalently, the same one-step
 *   transition of Frangioni and Gentile, 2006) plus the addition of the
 *   cost \f$ f_t( p ) = \alpha_t p^2 + \beta_t p + \gamma_t \f$.
 *   Transitions are arcs of length one between
 *   \f$ ( on^{\tau - 1} , t - 1 ) \f$ and \f$ ( on^\tau , t ) \f$.
 *
 * - The OFF portion of the state-space collapses Wuijts's \f$ M_{down} \f$
 *   counter into a *single* layer (see c_off_ready / c_off_any below).
 *
 * The shared machinery (data loading, the convex piecewise-quadratic value
 * function type and operations, and the whole spinning-reserve model
 * including the residual-ramp on->on transition sliding_min_corr) lives in
 * the base class ThermalUnitDPSolverBase, from which this class derives;
 * only the run-length DP structure, its state and its Solver interface are
 * here.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
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

#include "ThermalUnitDPSolverBase.h"

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

class ThermalUnitExtDPSolver : public ThermalUnitDPSolverBase
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// DIAGNOSTIC: DP-model cost of a given all-on power trajectory P (energy +
 /// folded reserve reward via corr = reserve_reward(min(A,B))). Lets a caller
 /// evaluate the MILP's trajectory under the DP's own reward model.
 double eval_allon_cost( const std::vector< double > & P ) const;

 /// DIAGNOSTIC: for each t, the DP's best cost-so-far to reach an on-state at
 /// power P[t] (min over surviving states of f_F[t](P[t]); +INF if no state
 /// covers P[t]). Compare against the cumulative cost of a trajectory to find
 /// the first t where the DP fails to propagate that trajectory.
 std::vector< double > ff_at_traj( const std::vector< double > & P ) const;

 /// DIAGNOSTIC: dump every surviving on-state at time t (tau, domain,
 /// value at p)
 void dump_states_at( Index t , double p ) const;

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

/*--------------------------------------------------------------------------*/
 /// returns the schedule the DP has found as a ThermalUnitBlockSolution
 /** Returns the schedule the dynamic programming has found as a
  * ThermalUnitBlockSolution [see ThermalUnitBlock.h], filled straight out of
  * the data structures of the Solver rather than by writing it in the
  * Variable of the ThermalUnitBlock and having it read back from there: no
  * abstract representation is therefore required to exist, and the
  * ThermalUnitBlock is not written into at all. See the analogous method of
  * ThermalUnitDPSolver for which parts of the solution are saved. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// the Solution is filled from the data of the DP, not from the Variable

 [[nodiscard]] bool is_get_Solution_physical( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// recovers the schedule the DP has found
 /** Recovers the schedule the dynamic programming has found: the active
  * power \p p, the commitment \p u, the primary and secondary spinning
  * reserve \p pr and \p sr, the reactive power \p q, and whether the unit
  * is \p built at all. It is what both get_var_solution() and
  * get_Solution() write, respectively into the Variable of the
  * ThermalUnitBlock and into the Solution. */

 virtual void recover_schedule( std::vector< double > & p ,
                                std::vector< double > & u ,
                                std::vector< double > & pr ,
                                std::vector< double > & sr ,
                                std::vector< double > & q ,
                                bool & built ) const;

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
  start = 0 ,     ///< nothing computed (or parameters have changed)
  loaded_OK = 1 , ///< parameters loaded from ThermalUnitBlock (reserved)
  dp_OK = 2 ,     ///< forward DP has been run; f_best_cost is valid
  sol_OK = 3      ///< P[] and U[] have been assembled by backtracking
 };

/*--------------------------------------------------------------------------*/
 /// summary of a single \f$ F^\tau_t \f$, kept alongside the PQFun itself
 /** Produced by run_DP() right after each \f$ F^\tau_t \f$ has been built.
  * It is used both to finalise the best cost at the end of the horizon and
  * to support backtracking without re-scanning the piecewise
  * representation:
  *
  * - min_val : the minimum value of \f$ F^\tau_t( p ) \f$ over the whole
  *             domain; TUEDPINF if the slot is infeasible (empty F)
  *
  * - argmin_p : a minimiser \f$ p^* \f$ of \f$ F^\tau_t \f$ on the same
  *              domain, the \f$ p^*_t \f$ of eq. (16) of Wuijts et al.
  *              (2021) */

 struct OnSlot {
  double min_val;
  double argmin_p;
 };

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// read all the parameters from the ThermalUnitBlock
 /** Loads the shared data via load_common_parameters() and then resets the
  * run-length solver's own output/pipeline state. virtual so a derived solver
  * (e.g. NuclearUnitExtDPSolver) can load its extra data on top. */
 virtual void load_parameters( void );

 /// process the queue of Modifications
 void process_modifications( void );

 /// process one Modification; returns true if a full reload is required
 /** virtual so that a derived solver can intercept its own Modifications. */
 virtual bool guts_of_process_modifications( const p_Mod mod );

/*--------------------------------------------------------------------------*/

 /// read the fixed status of the Variable of the ThermalUnitBlock
 /** Reads which commitment (and design) Variable are fixed, translating the
  * commitment fixings into the nxt_off / nxt_on tables that run_DP() uses
  * to kill the incompatible DP states; throws if any Variable that the DP
  * cannot honor (active power, reserves, start-up, shut-down, reactive) is
  * fixed. */
 void load_fixings( void );

/*--------------------------------------------------------------------------*/

 /// run the full forward DP (ON and OFF layers together)
 /** virtual so that a derived solver can replace it with a state-augmented
  * variant (e.g. NuclearUnitExtDPSolver, which adds the modulation
  * lockout). */
 virtual void run_DP( void );

 /// reconstruct the commitment/power schedule by backtracking the DP
 /** virtual for the same reason as run_DP(). */
 virtual void build_solution( void );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // NOTE: the data loaded from the ThermalUnitBlock (power/ramp/cost bounds,
 // reactive box, spinning-reserve factors and prices, design and fixing
 // data), the PieceQuad / PQFun value-function type and all its operations,
 // the whole reserve model (reserve_alloc / build_reserve_discount /
 // sliding_min_corr / reserve_corr_argmin / ...) and the PQ storage pool
 // live in the base class ThermalUnitDPSolverBase and are inherited.

 // ---- DP state ------------------------------------------------------- //

 char stage;                  ///< computation stage (see stage_value)
 bool f_solved;               ///< true iff f_best_cost encodes a feasible
                              ///< schedule (finite cost); set by run_DP()
 double f_best_cost;          ///< optimal objective value after run_DP()

 // -- ON-side value functions ------------------------------------------ //

 /// surviving F^tau_t functions, in sparse parallel-vector layout
 /** f_F[ t ] and f_tau[ t ] have the same length and are sorted by tau
  * ascending: f_F[ t ][ i ] is the piecewise convex quadratic value
  * function for the (t, f_tau[t][i]) DP state. Only *reachable* and
  * *relevant* (not RRF+-dominated) tau values are stored. */
 std::vector< std::vector< PQFun > > f_F;
 std::vector< std::vector< Index > > f_tau;

 /// per-slot summary (min_val, argmin_p) parallel to f_F[t] / f_tau[t]
 std::vector< std::vector< OnSlot > > f_on;

 // -- allocation pooling for run_DP() ---------------------------------- //
 // The per-step build buffers m_new_* are swapped into
 // f_F[t] / f_tau[t] / f_on[t] instead of freshly allocated; the PQFun
 // storage pool (m_pqpool) and sliding_min()'s scratch (m_raw) are
 // inherited from the base class.
 std::vector< PQFun >  m_new_F;
 std::vector< Index >  m_new_tau;
 std::vector< OnSlot > m_new_on;

 // -- OFF-side scalars ------------------------------------------------- //

 /// c_off_ready[ t ] : min cost of a schedule off at t AND off for at least
 /// min_down_time consecutive instants (legal to restart at t+1).
 std::vector< double > c_off_ready;

 /// c_off_any[ t ] : min cost of a schedule off at t, regardless of how long.
 /// Used only at the end of the horizon.
 std::vector< double > c_off_any;

 /// v_shutdown[ h ] : cost of reaching the "long shutdown arc" at the end of
 /// time h. Only well defined for h < time_horizon - 1.
 std::vector< double > v_shutdown;

 /// optimal (tau, p) that achieves v_shutdown[h]: needed to backtrack through
 /// the long shutdown arc.
 std::vector< Index  > v_shutdown_tau;
 std::vector< double > v_shutdown_p;

 /// "origin time" of c_off_ready[ t ]: the time h whose shutdown produced it
 /// via the long shutdown arc (-1 if +INF or inherited from the initial off
 /// trail).
 std::vector< int > f_ready_pred;

 /// analogous to f_ready_pred[ t ] but for c_off_any[ t ].
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
