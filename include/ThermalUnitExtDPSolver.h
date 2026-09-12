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
 * The problem is, over the horizon \f$ \mathcal{T} = \{ 0 , \ldots , n - 1
 * \} \f$,
 * \f[
 *  \min \Big\{ \, c( u ) + \sum_{t \in \mathcal{T}} f_t( p_t ) \, : \,
 *  P^{min}_t u_t \leq p_t \leq P^{max}_t u_t \, , \,
 *  p_t \leq p_{t-1} + u_{t-1} \Delta^+_{t-1} + ( 1 - u_{t-1} ) SU_t \, , \,
 *  p_{t-1} \leq p_t + u_t \Delta^-_{t-1} + ( 1 - u_t ) SD_t \, , \,
 *  u \in U \, \Big\}
 * \f]
 * where \f$ f_t( p ) = \alpha_t p^2 + \beta_t p + \gamma_t \f$ is the
 * convex production cost (\f$ \alpha_t \geq 0 \f$, while \f$ \beta_t \f$
 * may carry a Lagrangian price and be of any sign), \f$ c( u ) \f$ collects
 * the start-up costs, \f$ SU_t \f$ and \f$ SD_t \f$ are the start-up and
 * shut-down limits (bound_on and bound_down) and \f$ U \f$ is the set of
 * commitments that satisfy the minimum up- and down-time \f$ \tau^+ \f$ and
 * \f$ \tau^- \f$ from the initial state.
 *
 * <b>States.</b> For each \f$ t \f$, each run-length \f$ \tau \f$ (number of
 * consecutive time instants for which the unit has been on, ending with
 * \f$ t \f$) and each label \f$ \ell \f$, a *convex piecewise quadratic
 * function*
 * \f[ F^{\tau,\ell}_t : [ P^{min}_t , P^{max}_t ] \rightarrow \mathbb{R} \f]
 * stores the optimal cost of a schedule that is on at \f$ t \f$ with power
 * \f$ p \f$, whose current on-run has length exactly \f$ \tau \f$ and whose
 * state has label \f$ \ell \f$. The label is a small integer recording the
 * extra history a unit derived from ThermalUnitBlock needs for its own
 * temporal constraints: a thermal unit needs none and all its states have
 * label 0, a nuclear unit labels them by how long it is still locked out
 * from modulating [see NuclearUnitExtDPSolver]. The off-side collapses the
 * run-length of Wuijts et al. into two values per instant: \f$
 * c^{rdy}_t( e ) \f$, the optimal cost of a schedule off at \f$ t \f$ for at
 * least \f$ \tau^- \f$ instants (hence free to restart at \f$ t + 1 \f$)
 * whose off-state has label \f$ e \f$, and \f$ c^{any}_t \f$, the optimal
 * cost of a schedule off at \f$ t \f$ regardless of for how long, used only
 * at the end of the horizon.
 *
 * <b>Moves.</b> The on-states with label \f$ \ell \f$ at \f$ t - 1 \f$
 * continue at \f$ t \f$ by the moves \f$ m \in M_t( \ell ) \f$ of
 * on_moves(): each move has a landing label \f$ \ell_m \f$, a window
 * \f$ -w^-_m \leq p_t - p_{t-1} \leq w^+_m \f$ of the scheduled move, a
 * constant cost \f$ c_m \f$ and a range \f$ [ l_m , h_m ] \f$ of the landing
 * power. A thermal unit has the single move of window \f$ [ -\Delta^-_{t-1}
 * , \Delta^+_{t-1} ] \f$; the labels of the off-states evolve by
 * shut_label() \f$ \eta_t \f$ at a shut-down, idle_label() \f$ \iota_t \f$
 * along the idle instants and start_label() \f$ \sigma_t \f$ at a restart.
 *
 * <b>Recurrences.</b> With \f$ \widehat{f}_t = f_t + g^0_t \f$ the production
 * cost plus the capacity reward of the spinning reserve of a start-up
 * instant, whose band is capped by \f$ SU_t \f$ [see
 * ThermalUnitDPSolverBase::build_reserve_discount()],
 * \f[
 * \begin{array}{lll}
 *  \mbox{(on} \rightarrow \mbox{on)} &
 *  F^{\tau,\ell_m}_t( p ) = f_t( p ) + c_m + \min \{ F^{\tau-1,\ell}_{t-1}( q
 *  ) + corr_t( q , p ) \, : \, p - w^+_m \leq q \leq p + w^-_m \} &
 *  \tau > 1 \, , \, m \in M_t( \ell ) \, , \, p \in [ \max\{ P^{min}_t ,
 *  l_m \} , \min\{ P^{max}_t , h_m \} ] \\
 *  \mbox{(off} \rightarrow \mbox{on)} &
 *  F^{1,\sigma_t( e )}_t( p ) = \widehat{f}_t( p ) + SUC_t +
 *  c^{rdy}_{t-1}( e ) & p \in [ P^{min}_t , \min\{ P^{max}_t , SU_t \} ] \\
 *  \mbox{(on} \rightarrow \mbox{off)} &
 *  v^{sd}_h( e ) = \min \{ F^{\tau,\ell}_h( p ) \, : \, \tau \geq \tau^+ \, ,
 *  \, \eta_h( \ell ) = e \, , \, p \in [ P^{min}_h , SD_{h+1} ] \} &
 *  \mbox{shut-down at the end of } h \\
 *  \mbox{(off} \rightarrow \mbox{off)} &
 *  c^{rdy}_t( e ) = \min \{ c^{rdy}_{t-1}( e' ) \, : \, \iota_t( e' , 1 ) =
 *  e \} \cup \{ v^{sd}_{t-\tau^-}( e' ) \, : \, \iota_{t-\tau^-+1}( e' ,
 *  \tau^- ) = e \} &
 *  c^{any}_t = \min \{ c^{any}_{t-1} \, , \, \min_e v^{sd}_{t-1}( e ) \}
 * \end{array}
 * \f]
 * where the on->on step is the "ramp-constrained sliding minimum" of
 * Frangioni and Gentile (2006), taken over the window of the move, with the
 * penalty \f$ corr_t \f$ of the residual-ramp reserve folded in:
 * \f[
 *  corr_t( q , p ) = g_t\big( p , \min\{ A_t( p ) , B_t( p - q ) \} \big)
 *  \, , \quad A_t( p ) = \min\{ p - P^{min}_t , P^{max}_t - p \}
 *  \, , \quad B_t( d ) = \min\{ \Delta^+_{t-1} - d , \Delta^-_{t-1} + d \}
 * \f]
 * with \f$ g_t( p , H ) \leq 0 \f$ the greedy reward of the reserves held
 * in a band of width \f$ H \f$ [see ThermalUnitDPSolverBase::
 * sliding_min_corr()]. Note that the tent \f$ B_t \f$ of the reserve
 * deliverability is always the physical ramp of the step, even when the
 * window of the move is narrower. The long arc of the off->off step
 * enforces the minimum down-time in one jump, the initial off-trail of the
 * unit contributing a 0 at the instants in which it is already ready.
 * At a shut-down the reserve band of the closing instant is capped by
 * \f$ SD_{h+1} \f$ rather than \f$ P^{max}_h \f$, which requires re-running
 * its on->on step under that cap when a reserve is rewarded. The optimal
 * value is
 * \f[
 *  \min \big\{ \, c^{any}_{n-1} \, , \, \min_{\tau,\ell} \min_p
 *  F^{\tau,\ell}_{n-1}( p ) \, \big\}
 * \f]
 * plus the separable reactive term, and the schedule is recovered by
 * walking the recurrences backwards, each predecessor power being the
 * minimiser over the window of the move taken.
 *
 * <b>Pruning.</b> An \f$ F^{\tau,\ell}_t \f$ is dropped when some
 * \f$ F^{\tau',\ell'}_t \f$ with either \f$ \tau , \tau' \geq \tau^+ \f$ or
 * \f$ \tau = \tau' \f$, and with \f$ \ell' \f$ at least as good as \f$ \ell
 * \f$ (label_dominates()), is pointwise not larger: any schedule continuing
 * from the former is then matched, at no greater cost, by one continuing
 * from the latter, since the sliding minimum is monotone (Wuijts et al.,
 * Prop. 6.1). States with different labels are never merged into their
 * pointwise minimum, which would not be convex: each keeps its own
 * function. When trim_domination() the argument is used in full: the
 * domination is pointwise in \f$ p \f$, and a longer run-length reaches the
 * minimum up-time no later, so \f$ F^{\tau,\ell}_t \f$ loses the points of
 * its domain where some \f$ F^{\tau',\ell'}_t \f$ with \f$ \tau' \geq \tau
 * \f$ and \f$ \ell' \f$ at least as good is not larger (a union of
 * intervals, the difference of two quadratics changing sign at most twice),
 * what remains being one state per interval, each with the restriction of
 * the convex function. This is what keeps the states few when the domains
 * are narrow and shifted with respect to each other, as with a stable
 * output that does not change and moves at the full ramp. The surviving
 * states are few in practice, which makes the method run in close to
 * linear time, although its worst case is \f$ O( n^3 ) \f$ times the number
 * of labels.
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
 /// one on->on move out of an on-state, as returned by on_moves()
 /** The step from \f$ t - 1 \f$ to \f$ t \f$ of a unit that stays on:
  *
  * - lab : the label of the on-state the move lands in;
  *
  * - win_up, win_down : the window of the scheduled move,
  *   \f$ -win\_down \leq p_t - p_{t-1} \leq win\_up \f$; the tent of the
  *   reserve deliverability is always the physical ramp of the step, so a
  *   window narrower than the ramp (the modulation of a nuclear unit)
  *   limits the move without limiting the reserve;
  *
  * - cost : a constant cost of the move;
  *
  * - lo, hi : a range the landing power \f$ p_t \f$ is restricted to, on
  *   top of \f$ [ P^{min}_t , P^{max}_t ] \f$;
  *
  * - tag : an integer that the DP records, for each on instant of the
  *   optimal schedule, as the move taken to reach it [see U_move], and
  *   that a derived solver uses to recover the Variable it has on top of
  *   those of the ThermalUnitBlock.
  *
  * A window may exclude 0 (the move must then be strictly upwards or
  * downwards), and it may be degenerate, \f$ -win\_down = win\_up \f$: the
  * move is then exactly \f$ win\_up \f$, i.e., the value function is
  * shifted. */

 struct OnMove {
  Index lab;
  double win_up;
  double win_down;
  double cost;
  double lo;
  double hi;
  int tag;
 };

/*--------------------------------------------------------------------------*/
 /// how an on-state has been reached, for the backward pass
 /** Kept alongside each \f$ F^\tau_t \f$:
  *
  * - lab : the label of the on-state;
  *
  * - back : the index of the predecessor on-state in the (final) list at
  *   \f$ t - 1 \f$, BAD for a restart and for the states seeded at t = 0;
  *
  * - off : the label of the off-state the unit restarted from, for a
  *   restart (tau == 1);
  *
  * - move : the tag of the move that has been taken, -1 for a restart;
  *
  * - win_up, win_down : the window of that move. */

 struct OnLink {
  Index lab;
  std::size_t back;
  Index off;
  int move;
  double win_up;
  double win_down;
 };

 static constexpr std::size_t BAD = std::size_t( -1 );

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
  * fixed. virtual so that a derived solver can check its own Variable. */
 virtual void load_fixings( void );

/*--------------------------------------------------------------------------*/

 /// run the full forward DP (ON and OFF layers together)
 void run_DP( void );

 /// reconstruct the commitment/power schedule by backtracking the DP
 void build_solution( void );

/*--------------------------------------------------------------------------*/
/*------------------------- THE LABELS OF THE STATES -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name The labels of the states
 *
 * Each state of the DP carries, on top of the run-length \f$ \tau \f$ of
 * an on-state, a *label*: a small integer recording the extra history a
 * derived unit needs to enforce its own temporal constraints. A thermal
 * unit needs none, and all its states have label 0; a nuclear unit labels
 * its states by how long it is still locked out from modulating [see
 * NuclearUnitExtDPSolver]. The DP itself is the same for every unit: it
 * only asks the methods below which labels exist, which moves leave an
 * on-state, how a label evolves across the off-states and which label is
 * at least as good as another one. Every move is a sliding minimum over its
 * own window, so the value functions stay convex piecewise quadratic, and
 * the states with different labels are never merged: each keeps its own
 * function, and the domination pruning discards the redundant ones.
 *
 * The label of an on-state at \f$ t \f$ is the one *after* the decision
 * taken at \f$ t \f$, i.e., the one entering \f$ t + 1 \f$; the same holds
 * for the label of an off-state.
 *  @{ */

 /// the label that forbids a shut-down or a restart
 /** Returned by shut_label() when the unit cannot shut down after being on
  * with that label (e.g., in the middle of a modulation), and by
  * start_label() when it cannot restart with that off-label (e.g., when the
  * start-ups of the day are exhausted). */
 static constexpr Index NO_LABEL = Index( -1 );

 /// number of distinct labels of the on-states
 virtual Index on_labels( void ) const { return( 1 ); }

 /// number of distinct labels of the off-states
 virtual Index off_labels( void ) const { return( 1 ); }

 /// label of the state entering the time instant 0
 /** It is the label of the (on- or off-) state of the unit at the end of
  * the instant -1, the initial condition. */
 virtual Index init_label( void ) const { return( 0 ); }

 /// the moves out of an on-state with label \p lab at \p t - 1
 /** Appends to \p mv the moves the unit can take from an on-state with
  * label \p lab at \p t - 1 to an on-state at \p t (for \p t == 0 the
  * predecessor is the initial state, whose label is init_label()). The
  * default is the single move of a thermal unit, whose window is the ramp
  * of the step. */
 virtual void on_moves( Index t , Index lab ,
                        std::vector< OnMove > & mv ) const;

 /// label of the off-state of a unit shutting down after being on at \p t
 /** \p lab is the label of the on-state at \p t, the result the label of
  * the off-state at \p t + 1, or NO_LABEL if the shut-down is forbidden. */
 virtual Index shut_label( Index t , Index lab ) const { return( 0 ); }

 /// label of an off-state after \p k idle instants
 /** \p e is the label of the off-state entering \p t, the result the label
  * after the \p k idle instants \p t , ... , \p t + \p k - 1, i.e., the
  * one entering \p t + \p k. */
 virtual Index idle_label( Index t , Index e , Index k ) const {
  return( 0 );
  }

 /// label of the on-state of a unit restarting at \p t
 /** \p e is the label of the off-state entering \p t; NO_LABEL if the
  * restart is forbidden. */
 virtual Index start_label( Index t , Index e ) const { return( 0 ); }

/*--------------------------------------------------------------------------*/
 /// the labels of the on-states of a unit restarting at \p t
 /** Fills \p ls with the on-states that a unit restarting at \p t out of the
  * off-state with label \p e may land in: each entry is a label and the
  * range of the landing power that goes with it, and the ranges of two
  * entries only meet at their endpoints. The default is the single label of
  * start_label() over the whole range, which is what a unit whose rules do
  * not depend on the power it restarts at needs; a rule that does, such as
  * one whose labels say in which band of its range the output is, returns
  * one entry per band. */

 virtual void start_labels( Index t , Index e ,
                            std::vector< std::pair< Index ,
                                     std::pair< double , double > > > & ls )
  const {
  ls.clear();
  const Index lab = start_label( t , e );
  if( lab != NO_LABEL )
   ls.push_back( { lab , { - TUEDPINF , TUEDPINF } } );
  }

 /// true if label \p a is at least as good as label \p b for the future
 /** An on-state with label \p a can then take at least all the moves an
  * on-state with label \p b can, now and later, which is what allows the
  * former to prune the latter when its value function is pointwise not
  * larger. */
 virtual bool label_dominates( Index a , Index b ) const { return( true ); }

 /// true if the domination may also trim the domain of a state
 /** If true, an on-state loses the part of its domain where another one,
  * whose label and run-length allow it to prune it, is not larger (the
  * domination argument being pointwise in the power), which may split it
  * into two states. This matters when the domains of the value functions
  * are narrow and shifted with respect to each other, e.g., with moves at
  * the full ramp; for a thermal unit, whose value functions span the whole
  * reachable range, it is false. */
 virtual bool trim_domination( void ) const { return( false ); }

 /// the intervals of [ a , b ] where G is not larger than F
 static void not_larger( const PQFun & F , const PQFun & G , double a ,
                         double b ,
                         std::vector< std::pair< double , double > > & out );

/** @} ---------------------------------------------------------------------*/
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
 /** f_F[ t ], f_tau[ t ], f_on[ t ] and f_link[ t ] have the same length:
  * f_F[ t ][ i ] is the piecewise convex quadratic value function of the
  * on-state at t with run-length f_tau[ t ][ i ] and label
  * f_link[ t ][ i ].lab. Only *reachable* and *relevant* (not dominated)
  * states are stored. */
 std::vector< std::vector< PQFun > > f_F;
 std::vector< std::vector< Index > > f_tau;

 /// per-slot summary (min_val, argmin_p) parallel to f_F[t] / f_tau[t]
 std::vector< std::vector< OnSlot > > f_on;

 /// per-slot label and predecessor, parallel to f_F[t] / f_tau[t]
 std::vector< std::vector< OnLink > > f_link;

 // -- allocation pooling for run_DP() ---------------------------------- //
 // The per-step build buffers m_new_* are swapped into
 // f_F[t] / f_tau[t] / f_on[t] / f_link[t] instead of freshly allocated;
 // the PQFun storage pool (m_pqpool) and sliding_min()'s scratch (m_raw)
 // are inherited from the base class.
 std::vector< PQFun >  m_new_F;
 std::vector< Index >  m_new_tau;
 std::vector< OnSlot > m_new_on;
 std::vector< OnLink > m_new_link;
 std::vector< OnMove > m_moves;

 /// the labels a restart may land in, with their ranges [see start_labels()]
 mutable std::vector< std::pair< Index , std::pair< double , double > > >
  m_start_labs;   ///< scratch for on_moves()
 /// scratch of the domination with trimming
 std::vector< std::pair< double , double > > m_cover , m_rest;

 // -- OFF-side values -------------------------------------------------- //
 // the arrays indexed by time and off label are flat, entry t * E + e with
 // E = off_labels(), so that a unit with a single label keeps one value per
 // time instant

 /// c_off_ready[ t * E + e ] : min cost of a schedule off at t AND off for at
 /// least min_down_time consecutive instants (legal to restart at t+1), with
 /// label e entering t+1.
 std::vector< double > c_off_ready;

 /// c_off_any[ t ] : min cost of a schedule off at t, regardless of how long
 /// and of the label. Used only at the end of the horizon.
 std::vector< double > c_off_any;

 /// v_shutdown[ h * E + e ] : cost of reaching the "long shutdown arc" at the
 /// end of time h, with label e entering h+1. Only well defined for
 /// h < time_horizon - 1.
 std::vector< double > v_shutdown;

 /// optimal (tau, p, link) that achieves v_shutdown[ h * E + e ]: needed to
 /// backtrack through the long shutdown arc. The link is that of the
 /// closing on-state at h, which may not be in f_F[ h ] when the closing
 /// transition has been re-run under the shut-down cap.
 std::vector< Index  > v_shutdown_tau;
 std::vector< double > v_shutdown_p;
 std::vector< OnLink > v_shutdown_link;

 /// "origin time" of c_off_ready[ t * E + e ]: the time h whose shutdown
 /// produced it via the long shutdown arc (-1 if +INF or inherited from the
 /// initial off trail), and the label entering h+1 of that shutdown.
 std::vector< int > f_ready_pred;
 std::vector< Index > f_ready_lab;

 /// analogous to f_ready_pred / f_ready_lab but for c_off_any[ t ].
 std::vector< int > f_any_pred;
 std::vector< Index > f_any_lab;

 // ---- output ---------------------------------------------------------- //

 /// power profile of the optimal schedule, as filled by build_solution()
 std::vector< double > P;
 /// commitment profile of the optimal schedule, as filled by build_solution()
 std::vector< bool >   U;
 /// for each on instant of the optimal schedule, the label of the on-state
 /// and the tag of the move taken to reach it (-1 for a restart), as filled
 /// by build_solution(); meaningless where the unit is off
 std::vector< Index > U_lab;
 std::vector< int >   U_move;

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( ThermalUnitExtDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitExtDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File ThermalUnitExtDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
