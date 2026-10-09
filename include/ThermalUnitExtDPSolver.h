/*--------------------------------------------------------------------------*/
/*---------------------- File ThermalUnitExtDPSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitExtDPSolver class, a Solver of the
 * single-unit commitment problem of a ThermalUnitBlock by a dynamic program
 * whose on-states are indexed by the length of the current on-run, while the
 * off-states are collapsed into one value per instant (and per label).
 *
 * We solve the problem of ThermalUnitDPSolverBase, i.e., that of the rows of
 * ThermalUnitBlock::generate_abstract_constraints(). In the file comment of
 * ThermalUnitDPSolverBase one finds the model, the notation, the cost
 * \f$ f_t( p ) = a_t p^2 + b_t p + c_t \f$ of an on instant and the reserve
 * term \f$ g_t( p , \mathcal{H}_t( p , q ) ) \f$ of (2)-(3) there.
 *
 * <b>States.</b> For each instant \f$ t \f$, run-length \f$ \tau \geq 1 \f$
 * (the number of consecutive instants the unit has been on, \f$ t \f$
 * included) and label \f$ \ell \f$, the convex piecewise quadratic function
 * \f$ F^{\tau,\ell}_t( p ) \f$ is the least cost of a schedule of the
 * instants \f$ 0 , \ldots , t \f$ that is on at \f$ t \f$ with
 * \f$ p^{ac}_t = p \f$, whose current on-run has length \f$ \tau \f$ and
 * whose state after the decision at \f$ t \f$ has label \f$ \ell \f$. The
 * label is a small integer that records the history which a unit derived from
 * ThermalUnitBlock needs for its rules. For instance, a thermal unit has the
 * single label 0, while a nuclear unit labels its states by its modulation
 * lockout, its mode, the band of its output and the counters of the day [see
 * NuclearUnitExtDPSolver]. On the off side there are two values per instant.
 * The first is \f$ c^{rdy}_t( o ) \f$, the least cost of a schedule that is
 * off at \f$ t \f$ since at least \f$ \tau^- \f$ instants (hence free to
 * start up at \f$ t + 1 \f$) and has off-label \f$ o \f$; the second is
 * \f$ c^{any}_t \f$, the least cost of a schedule that is off at \f$ t \f$,
 * which is used only at the end of the horizon. Finally, the values
 * \f$ \phi^{sd}_h( o ) \f$ are the least costs of the schedules that are on
 * at \f$ h \f$ and off at \f$ h + 1 \f$ with off-label \f$ o \f$.
 *
 * <b>Moves.</b> An on-state with label \f$ \ell \f$ at \f$ t - 1 \f$
 * continues at \f$ t \f$ by the moves \f$ j \in J_t( \ell ) \f$ of
 * on_moves(): each has a landing label \f$ \ell_j \f$, a window
 * \f$ -\omega^-_j \leq p^{ac}_t - p^{ac}_{t-1} \leq \omega^+_j \f$ of the
 * scheduled move, a constant cost \f$ c_j \f$ and a range
 * \f$ [ p^{lo}_j , p^{hi}_j ] \f$ of the power at \f$ t \f$. A thermal unit
 * has the single move with \f$ \omega^+_j = \Delta^+_t \f$,
 * \f$ \omega^-_j = \Delta^-_t \f$ (the ramps of the step from \f$ t - 1 \f$
 * to \f$ t \f$), \f$ c_j = 0 \f$ and no range. In turn, the off-labels evolve
 * by shut_label() \f$ \varsigma_h( \ell ) \f$ at a shut-down after \f$ h \f$
 * and by idle_label() \f$ \iota_t( o , k ) \f$ along \f$ k \f$ idle instants
 * from \f$ t \f$. A start-up at \f$ t \f$ from the off-label \f$ o \f$ may
 * land in the on-labels of the set \f$ S_t( o ) \f$ given by start_labels(),
 * each with its range of the power. Also, a label NO_LABEL forbids the
 * shut-down or the start-up.
 *
 * <b>Recursion.</b> For \f$ t \geq 1 \f$,
 * \f{align*}{
 *   F^{\tau,\ell_j}_t( p ) &= f_t( p ) + c_j + \min \bigl\{ \,
 *     F^{\tau-1,\ell}_{t-1}( q ) + g_t( p , \mathcal{H}_t( p , q ) ) \, : \,
 *     p - \omega^+_j \leq q \leq p + \omega^-_j \, \bigr\} \, , \quad
 *     \tau \geq 2 \, , \; j \in J_t( \ell ) \, , \tag{1} \\
 *   F^{1,\ell'}_t( p ) &= f_t( p ) + c^{su}_t + g_t( p , \min\{ p -
 *     P^{mn}_t , P^{su}_t - p \} ) + \min \{ \, c^{rdy}_{t-1}( o ) \, : \,
 *     \ell' \in S_t( o ) \, \} \, , \tag{2} \\
 *   \phi^{sd}_h( o ) &= c^{sd}_{h+1} + \min \{ \, \widetilde F^{\tau,\ell}_h(
 *     p ) \, : \, \tau \geq \tau^+ \, , \; \varsigma_h( \ell ) = o \, , \;
 *     P^{mn}_h \leq p \leq P^{sd}_{h+1} \, \} \, , \tag{3} \\
 *   c^{rdy}_t( o ) &= \min \Bigl( \bigl\{ \, c^{rdy}_{t-1}( o' ) \, : \,
 *     \iota_t( o' , 1 ) = o \, \bigr\} \cup \bigl\{ \,
 *     \phi^{sd}_{t-\tau^-}( o' ) \, : \, \iota_{t-\tau^-+1}( o' , \tau^- ) = o
 *     \, \bigr\} \Bigr) \, , \tag{4} \\
 *   c^{any}_t &= \min \bigl\{ \, c^{any}_{t-1} \, , \; \min_o
 *     \phi^{sd}_{t-1}( o ) \, \bigr\} \, , \tag{5}
 * \f}
 * where (1) is defined for \f$ p \in [ \max\{ P^{mn}_t , p^{lo}_j \} , \min\{
 * P^{mx}_t , p^{hi}_j \} ] \f$, (2) for \f$ p \in [ P^{mn}_t , \min\{
 * P^{mx}_t , P^{su}_t \} ] \f$ intersected with the range of \f$ \ell' \f$
 * (the start-up cost \f$ c^{su}_t \f$ depends only on the instant of the
 * start-up, since the state does not carry how long the unit has been off).
 * In (3), \f$ \widetilde F^{\tau,\ell}_h \f$ is \f$ F^{\tau,\ell}_h \f$ with
 * the cap of the reserve band of instant \f$ h \f$ lowered to
 * \f$ P^{sd}_{h+1} \f$; this requires recomputing (1) at \f$ h \f$ (or
 * adjusting (2), for \f$ \tau = 1 \f$) when a reserve cost is negative and
 * the cap is below \f$ P^{mx}_h \f$. The minimum in (1) is the transition (4)
 * of ThermalUnitDPSolverBase, computed by
 * ThermalUnitDPSolverBase::sliding_min_corr() (exactly, save for the parts
 * that it interpolates; see there). Note that the ramp terms of the reserve
 * band are always those of the ramp of the step, also when the window of the
 * move is narrower. The second set in (4) is the shut-down at
 * \f$ h = t - \tau^- \f$ followed by \f$ \tau^- \f$ idle instants, which
 * enforces the minimum down time in one jump. In (3), the cost
 * \f$ c^{sd}_{h+1} \f$ is 0 for \f$ h = T - 1 \f$, and
 * \f$ \phi^{sd}_{T-1} \f$ is not used. A state at \f$ t \f$ does not exist if
 * the commitment Variable is fixed to the opposite value at \f$ t \f$, and
 * the second set in (4) needs no instant fixed on in
 * \f$ h + 1 , \ldots , t \f$ [see load_fixings()].
 *
 * <b>Initial state.</b> If the unit is on before the horizon
 * (\f$ \tau_0 > 0 \f$), the on-states at \f$ t = 0 \f$ are
 * \f$ F^{\tau_0+1,\ell_j}_0( p ) = f_0( p ) + c_j + g_0( p , \mathcal{H}_0( p
 * , p_{-1} ) ) \f$ for the moves \f$ j \f$ out of the initial label
 * init_label(), on \f$ [ \max\{ P^{mn}_0 , p_{-1} - \omega^-_j , p^{lo}_j \}
 * , \min\{ P^{mx}_0 , p_{-1} + \omega^+_j , p^{hi}_j \} ] \f$. Hence, the
 * reserve of instant 0 is bounded by the ramp left after the move from
 * InitialPower, as in the deliverability rows of ThermalUnitBlock at 0 [see
 * ThermalUnitDPSolverBase::initial_reserve_discount()]. Moreover, if
 * \f$ \tau_0 \geq \tau^+ \f$, the unit may also be off at 0, at the cost
 * \f$ c^{sd}_0 \f$ of the shut-down at 0, unless either
 * \f$ p_{-1} > P^{sd}_0 \f$, which all the formulations of ThermalUnitBlock
 * forbid whether or not DeltaRampDown is given (see
 * ThermalUnitBlock::generate_abstract_constraints()), or DeltaRampUp is given
 * and \f$ p_{-1} < P^{mn}_0 \f$, which ThermalUnitBlock excludes. When it is
 * off there, the unit is then ready to start up at \f$ t + 1 \f$ as soon as
 * \f$ t + 1 \geq \tau^- \f$, which adds the value \f$ c^{sd}_0 \f$ (the
 * initial off-trail) to the sets of (4). If, instead, the unit is off before
 * the horizon (\f$ \tau_0 \leq 0 \f$), it may start up at 0 if
 * \f$ - \tau_0 \geq \tau^- \f$, by (2) with \f$ c^{rdy}_{-1} = 0 \f$. It may
 * also be off at 0 at no cost, and it is then ready to start up at
 * \f$ t + 1 \f$ as soon as \f$ - \tau_0 + t + 1 \geq \tau^- \f$, which adds
 * the value 0 to the sets of (4). In both cases the initial off-trail needs
 * no instant fixed on in \f$ 0 , \ldots , t \f$.
 *
 * <b>End of the horizon.</b> The optimal value is
 * \f[
 *   \min \Bigl\{ \, c^{any}_{T-1} \, , \; \min_{ \tau , \ell } \min_p
 *     F^{\tau,\ell}_{T-1}( p ) \, \Bigr\} \tag{6}
 * \f]
 * (no shut-down limit applies at \f$ T - 1 \f$, since the unit is not shut
 * down within the horizon), plus the constant of the reactive power when the
 * unit has it. If the unit has an InvestmentCost, it is built if this value
 * plus the current coefficient of the design variable is not positive (or if
 * a fixing forces it), and the value is 0 otherwise. Also, the reported value
 * is multiplied by the scale factor \f$ \sigma \f$. Finally, the schedule is
 * recovered by walking the recursion backwards, taking as power at
 * \f$ t - 1 \f$ a minimizer in (1) at the power at \f$ t \f$ [see
 * ThermalUnitDPSolverBase::reserve_corr_argmin()], and the reserves are those
 * of (2) of ThermalUnitDPSolverBase at the recovered schedule.
 *
 * <b>Pruning.</b> A state \f$ ( \tau , \ell , F ) \f$ at \f$ t \f$ is
 * discarded if another state \f$ ( \tau' , \ell' , F' ) \f$ at \f$ t \f$ has
 * \f$ F' \leq F \f$ everywhere on the domain of \f$ F \f$, a label
 * \f$ \ell' \f$ at least as good as \f$ \ell \f$ [see label_dominates()], and
 * either \f$ \tau , \tau' \geq \tau^+ \f$ or \f$ \tau' = \tau \f$ (or, if
 * trim_domination(), \f$ \tau' \geq \tau \f$). This rule is exact. In fact,
 * the moves, the start-ups and the shut-downs out of a state depend only on
 * \f$ t \f$, on its label and on the power, and on the run-length only
 * through the condition \f$ \tau \geq \tau^+ \f$ of (3), which \f$ \tau' \f$
 * satisfies whenever \f$ \tau \f$ does. Furthermore, a label at least as good
 * allows any sequence of moves of the other one, landing in labels that are
 * again at least as good, and (1) is monotone in
 * \f$ F^{\tau-1,\ell}_{t-1} \f$. Hence, any completion of a schedule through
 * \f$ ( \tau , \ell , p ) \f$ is also a completion through
 * \f$ ( \tau' , \ell' , p ) \f$, at no greater cost, and the optimal value
 * (6) does not change. With trim_domination() the argument (which is
 * pointwise in \f$ p \f$) is also applied to the part of the domain of
 * \f$ F \f$ where \f$ F' \f$ is not larger, which is a union of intervals
 * since the difference of two quadratic pieces changes sign at most twice.
 * That part is removed, and what remains becomes one state per interval, each
 * with the restriction of \f$ F \f$; this keeps the states few when the
 * domains are narrow and shifted with respect to each other. Note that states
 * with different labels are never merged into their pointwise minimum, which
 * would not be convex. The surviving states are typically few, and therefore
 * the method often runs in a time close to linear in \f$ T \f$; however, in
 * the worst case their number grows linearly with \f$ t \f$ (times the number
 * of labels).
 *
 * The minimum up and down times are taken to be at least 1, as in
 * ThermalUnitBlock, and therefore the value (6) is the optimal value of the
 * rows of ThermalUnitBlock, with the exceptions stated in
 * ThermalUnitDPSolverBase.h (a ReferenceSchedule and the fixed Variable other
 * than the commitment and the design, which are refused, and the parts of the
 * transitions computed by interpolation).
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
/// dynamic programming Solver of a ThermalUnitBlock over the run-lengths

class ThermalUnitExtDPSolver : public ThermalUnitDPSolverBase
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// diagnostic: the cost, in the model of the DP, of the schedule on at all
 /// the instants with active power P, the reserve term included
 double eval_allon_cost( const std::vector< double > & P ) const;

 /// diagnostic: for each t, the least value at P[ t ] of the on-states that
 /// survive at t, +INF if none covers P[ t ]
 std::vector< double > ff_at_traj( const std::vector< double > & P ) const;

 /// diagnostic: print the on-states that survive at t, with their
 /// run-length, domain and value at p
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
 OFValue get_lb( void ) override { return( scaled_value( f_best_cost ) ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( scaled_value( f_best_cost ) ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override {
  return( scaled_value( f_best_cost ) );
  }

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
 /// the minimum of an on-state, kept next to its value function
 /** For an on-state \f$ F^{\tau,\ell}_t \f$, as built by run_DP():
  *
  * - min_val: the minimum of \f$ F^{\tau,\ell}_t \f$ over its domain,
  *   TUEDPINF if the state is infeasible (empty function);
  *
  * - argmin_p: a minimizer of \f$ F^{\tau,\ell}_t \f$, the power at
  *   \f$ t \f$ the backward pass starts from if the state ends the optimal
  *   schedule. */

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
  * - win_up, win_down : the window \f$ \omega^+_j \f$, \f$ \omega^-_j \f$ of
  *   the scheduled move, \f$ -\omega^-_j \leq p^{ac}_t - p^{ac}_{t-1} \leq
  *   \omega^+_j \f$; the ramp terms of the reserve band are always those of
  *   the ramp of the step, so that a window narrower than the ramp limits the
  *   move without limiting the reserve;
  *
  * - cost : a constant cost of the move;
  *
  * - lo, hi : the range \f$ [ p^{lo}_j , p^{hi}_j ] \f$ the power \f$ p^{ac}_t
  *   \f$ is restricted to, on top of \f$ [ P^{mn}_t , P^{mx}_t ] \f$;
  *
  * - tag : an integer that the DP records, for each on instant of the
  *   optimal schedule, as the move taken to reach it [see U_move], and
  *   that a derived solver uses to recover the Variable it has on top of
  *   those of the ThermalUnitBlock.
  *
  * A window may exclude 0 (the move must then go upwards, or downwards),
  * and it may be a single point, \f$ -\omega^-_j = \omega^+_j \f$: the move is
  * then
  * exactly \f$ \omega^+_j \f$, i.e., the value function is translated. */

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
 /** Kept next to each \f$ F^{\tau,\ell}_t \f$:
  *
  * - lab : the label of the on-state;
  *
  * - back : the index of the predecessor on-state in the list at
  *   \f$ t - 1 \f$, BAD for a start-up and for the states at
  *   \f$ t = 0 \f$;
  *
  * - off : the label of the off-state the unit started up from, for a
  *   start-up (\f$ \tau = 1 \f$);
  *
  * - move : the tag of the move that has been taken, -1 for a start-up;
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
  * state of this solver; virtual, so that a derived solver (e.g.,
  * NuclearUnitExtDPSolver) can load its own data on top. */
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
  * to kill the incompatible DP states; throws std::logic_error if any
  * Variable that the DP cannot honor (active power, reserves, start-up,
  * shut-down, reactive, or a Variable of an extended formulation [see
  * fixed_extended_variable()]) is fixed. virtual so that a derived solver
  * can check its own Variable [see reads_group()]. */
 virtual void load_fixings( void );

/*--------------------------------------------------------------------------*/

 /// run the forward recursion (1)-(5) of the file comment
 void run_DP( void );

 /// recover the optimal schedule by walking the recursion backwards
 void build_solution( void );

/*--------------------------------------------------------------------------*/
/*------------------------- THE LABELS OF THE STATES -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name The labels of the states
 *
 * Each state of the DP carries, on top of the run-length \f$ \tau \f$ of
 * an on-state, a label: a small integer recording the extra history a
 * derived unit needs to enforce its own temporal constraints. A thermal
 * unit needs none, and all its states have label 0; a nuclear unit labels
 * its states by its mode, by how long it is still locked out from
 * modulating, by the band of its output and by the counters of the day
 * [see NuclearUnitExtDPSolver]. The DP itself is the same for every unit: it
 * only asks the methods below which labels exist, which moves leave an
 * on-state, how a label evolves across the off-states and which label is
 * at least as good as another one. Every move is a sliding minimum over its
 * own window, so the value functions stay convex piecewise quadratic, and
 * the states with different labels are never merged: each keeps its own
 * function, and the domination pruning discards the redundant ones.
 *
 * The label of an on-state at \f$ t \f$ is the one after the decision
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
  * \f$ [ -\Delta^-_t , \Delta^+_t ] \f$ of the step. */
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
 /** An on-state with label \p a can then take at least all the moves, the
  * start-ups after a shut-down and the shut-downs an on-state with label
  * \p b can, now and later, landing in labels that are again at least as
  * good, which is what allows the former to prune the latter when its value
  * function is pointwise not larger [see the file comment]. */
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

 /// the on-states that survive at each instant
 /** f_F[ t ], f_tau[ t ], f_on[ t ] and f_link[ t ] have the same length:
  * f_F[ t ][ i ] is the value function of the on-state at t with run-length
  * f_tau[ t ][ i ] and label f_link[ t ][ i ].lab. Only the states that are
  * reachable and that the pruning keeps are stored. */
 std::vector< std::vector< PQFun > > f_F;
 std::vector< std::vector< Index > > f_tau;

 /// the minimum of each on-state, parallel to f_F[ t ]
 std::vector< std::vector< OnSlot > > f_on;

 /// the label and the predecessor of each on-state, parallel to f_F[ t ]
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

 /// c_off_ready[ t * E + e ] is \f$ c^{rdy}_t( o ) \f$, o = e, of the file
 /// comment
 std::vector< double > c_off_ready;

 /// c_off_any[ t ] is \f$ c^{any}_t \f$ of the file comment
 std::vector< double > c_off_any;

 /// v_shutdown[ h * E + e ] is \f$ \phi^{sd}_h( o ) \f$, o = e, of the file
 /// comment, defined for h < time_horizon - 1
 std::vector< double > v_shutdown;

 /// the run-length, the power and the link of the on-state at h that
 /// attains v_shutdown[ h * E + e ], for the backward pass; the link is that
 /// of the closing on-state at h, which is not in f_F[ h ] when (1) has been
 /// recomputed at h under the shut-down cap
 std::vector< Index  > v_shutdown_tau;
 std::vector< double > v_shutdown_p;
 std::vector< OnLink > v_shutdown_link;

 /// the instant h of the shut-down that c_off_ready[ t * E + e ] comes from
 /// through the second set of (4), -1 if it is +INF or it comes from the
 /// initial off-trail, and the off-label of that shut-down
 std::vector< int > f_ready_pred;
 std::vector< Index > f_ready_lab;

 /// the same as f_ready_pred and f_ready_lab, for c_off_any[ t ]
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
