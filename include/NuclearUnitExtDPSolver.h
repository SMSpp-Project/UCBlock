/*--------------------------------------------------------------------------*/
/*--------------------- File NuclearUnitExtDPSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the NuclearUnitExtDPSolver class, a Solver of the
 * single-unit commitment problem of a NuclearUnitBlock (i.e., of a thermal
 * unit subject to the rules on the modulation of its output) by the dynamic
 * program of ThermalUnitExtDPSolver, whose states carry a label that enforces
 * those rules.
 *
 * <b>The model.</b> On top of the rows of ThermalUnitBlock, a
 * NuclearUnitBlock has the binary variables \f$ m_t \f$ (a modulation step at
 * \f$ t \f$), \f$ d_t \f$ (the step goes downwards), \f$ \delta_t \f$ (a deep
 * decrease at \f$ t \f$) and \f$ b^k_t \f$ (the output at \f$ t \f$ is in
 * band \f$ k \f$), and the rows of
 * NuclearUnitBlock::generate_abstract_constraints() and of its operating
 * rules. In the notation of NuclearUnitBlock, these rows are the following:
 *
 * - An on instant \f$ t \f$ that is not a start-up is either stable or a
 *   modulation step. In the former case the output moves within
 *   \f$ [ -\Delta^{M-}_t , \Delta^{M+}_t ] \f$ (the ramps
 *   ModulationDeltaRampDown and ModulationDeltaRampUp, at most the ramps
 *   \f$ \Delta^-_t \f$ and \f$ \Delta^+_t \f$ of the step); in the latter
 *   \f$ m_t = 1 \f$, and the output moves within the ramps
 *   \f$ [ -\Delta^-_t , \Delta^+_t ] \f$, downwards if \f$ d_t = 1 \f$ and
 *   upwards otherwise when the direction matters.
 *
 * - A modulation is a sequence of at most \f$ L^M \f$ consecutive steps
 *   (MaxModulationLength) in the same direction, in which all the steps but
 *   the last move the output by exactly the ramp, and no step follows its
 *   last step for \f$ \tau^M - 1 \f$ instants (ModulationTime
 *   \f$ \tau^M \geq 2 \f$). There is no step at
 *   \f$ t < \tau^M - \tau^M_0 \f$, where the unit has modulated last
 *   \f$ \tau^M_0 \geq 1 \f$ instants before 0 (InitModulation), nor at the
 *   \f$ \tau^v \f$ instants that begin with a start-up
 *   (StabilityAfterStartUp).
 *
 * - In each day, i.e., in each of the disjoint sets of \f$ T^{day} \f$
 *   consecutive instants (DayLength, the whole horizon if 0), at most
 *   \f$ N^M \f$ modulations start (ModulationsPerDay), at most \f$ N^{dd} \f$
 *   deep decreases happen (DeepDecreasesPerDay) and at most \f$ N^{su} \f$
 *   start-ups happen (StartUpsPerDay). Each limit applies only if it is
 *   given.
 *
 * - A deep decrease happens at an on instant \f$ t \f$ with the unit on at
 *   \f$ t - 1 \f$ when \f$ p^{ac}_{t-1} - p^{ac}_t > \tilde\Delta_t \f$ and
 *   \f$ p^{ac}_t < \tilde p_t \f$ (DeepDecreaseGradient,
 *   DeepDecreaseThreshold), at the cost \f$ c^{dd}_t \geq 0 \f$; each
 *   downward step costs \f$ c^-_t \geq 0 \f$ (DownModulationCost).
 *
 * - If the output is split into three bands by the breakpoints
 *   \f$ P^b_1 < P^b_2 \f$ (PowerBands), the output of a stable instant stays
 *   in the band of the previous one, while a modulation crosses exactly one
 *   boundary of its band, towards the adjacent band in its direction. In
 *   particular, the steps that precede the last one keep the output in the
 *   band of origin, and its last step lands in the adjacent band. Hence, with
 *   a constant full ramp \f$ \Delta \f$, a modulation of \f$ k \f$ steps that
 *   starts from the output \f$ p \f$ and crosses the breakpoint \f$ P^b \f$
 *   needs \f$ ( k - 1 ) \Delta \leq | p - P^b | \leq k \Delta \f$; therefore,
 *   it lasts \f$ \lceil | p - P^b | / \Delta \rceil \f$ steps, or that number
 *   or one more if the ratio is an integer (the last step of the longer one
 *   moves the output by 0 from \f$ P^b \f$; see "PowerBands" in
 *   NuclearUnitBlock::deserialize()), and at most \f$ L^M \f$. Also, no
 *   modulation goes upwards from the highest band nor downwards from the
 *   lowest one, and the band at \f$ -1 \f$ is the lowest one that contains
 *   InitialPower.
 *
 * - The unit may shut down at any instant but in the middle of a modulation,
 *   i.e., after a step that is not the last one; in particular, it may shut
 *   down during the stability that follows a modulation or a start-up.
 *
 * - Also in the bands, the bounds of the output are the operational ones of
 *   ThermalUnitBlock; hence, at an instant with availability 0 the unit may
 *   be on at the output 0, in the lowest band, as in the rows of
 *   NuclearUnitBlock.
 *
 * <b>The labels.</b> The label of an on-state entering \f$ t + 1 \f$ is
 * the tuple \f$ ( \varpi , n^{lk} , k , n^M , n^{dd} , n^{su} ) \f$. Its
 * first component is the mode \f$ \varpi \f$, which is \f$ \mathrm{st} \f$
 * (stable, i.e., not in the middle of a modulation), \f$ \mathrm{up} \f$ or
 * \f$ \mathrm{dn} \f$ (in the middle of an upward or downward modulation).
 * The second one is, for \f$ \varpi = \mathrm{st} \f$, the lockout \f$ n^{lk}
 * \in \{ 0 , \ldots , \max\{ \tau^M , \tau^v \} - 1 \} \f$, i.e., the number
 * of instants from \f$ t + 1 \f$ on in which a step is still forbidden, and
 * for \f$ \varpi \in \{ \mathrm{up} , \mathrm{dn} \} \f$ the number
 * \f$ n^{lk} \in \{ 1 , \ldots , L^M - 1 \} \f$ of steps taken so far. Then
 * comes the band \f$ k \f$ of the output (for \f$ \varpi \in \{ \mathrm{up} ,
 * \mathrm{dn} \} \f$, the band the modulation left), which exists only if the
 * output is banded, and finally the numbers of modulations, deep decreases
 * and start-ups in the day of \f$ t + 1 \f$, each of which exists only if the
 * corresponding limit is given. For an off-state, the label is
 * \f$ ( n^{lk} , n^M , n^{dd} , n^{su} ) \f$. Let
 * \f$ n^{lk-} = \max\{ n^{lk} - 1 , 0 \} \f$, let \f$ \Pi_k \f$ be the range
 * of the output in band \f$ k \f$ (the closed interval between its
 * breakpoints, or everything if the output is not banded), and read
 * \f$ k \pm 1 \f$ as \f$ k \f$ when the output is not banded. The moves out
 * of an on-state with label \f$ ( \varpi , n^{lk} , k ) \f$ at \f$ t - 1 \f$
 * are then (omitting the counters where they do not change)
 * \f[
 *   \begin{array}{llllll}
 *     \mbox{from} & p^{ac}_t - p^{ac}_{t-1} & \mbox{to} & p^{ac}_t \in &
 *       \mbox{cost} & \mbox{when} \\
 *     ( \mathrm{st} , n^{lk} , k ) & [ -\Delta^{M-}_t , \Delta^{M+}_t ] &
 *       ( \mathrm{st} , n^{lk-} , k ) & \Pi_k & 0 & \\
 *     ( \mathrm{st} , 0 , k ) & [ -\Delta^-_t , \Delta^+_t ] &
 *       ( \mathrm{st} , \tau^M - 1 , k ) , \; n^M + 1 & & 0 & \mbox{(i)} \\
 *     ( \mathrm{st} , 0 , k ) & [ 0 , \Delta^+_t ] &
 *       ( \mathrm{st} , \tau^M - 1 , k + 1 ) , \; n^M + 1 & \Pi_{k+1} & 0 &
 *       \mbox{(ii)} \\
 *     ( \mathrm{st} , 0 , k ) & [ -\Delta^-_t , 0 ] &
 *       ( \mathrm{st} , \tau^M - 1 , k - 1 ) , \; n^M + 1 & \Pi_{k-1} &
 *       c^-_t & \mbox{(iii)} \\
 *     ( \mathrm{st} , 0 , k ) & \{ \Delta^+_t \} & ( \mathrm{up} , 1 , k ) ,
 *       \; n^M + 1 & \Pi_k & 0 & \mbox{(ii), (iv)} \\
 *     ( \mathrm{st} , 0 , k ) & \{ -\Delta^-_t \} & ( \mathrm{dn} , 1 , k ) ,
 *       \; n^M + 1 & \Pi_k & c^-_t & \mbox{(iii), (iv)} \\
 *     ( \mathrm{up} , n^{lk} , k ) & \{ \Delta^+_t \} & ( \mathrm{up} , n^{lk}
 *       + 1 , k ) & \Pi_k & 0 & n^{lk} + 1 < L^M \\
 *     ( \mathrm{up} , n^{lk} , k ) & [ 0 , \Delta^+_t ] &
 *       ( \mathrm{st} , \tau^M - 1 , k + 1 ) & \Pi_{k+1} & 0 & \\
 *     ( \mathrm{dn} , n^{lk} , k ) & \{ -\Delta^-_t \} & ( \mathrm{dn} ,
 *       n^{lk} + 1 , k ) & \Pi_k & c^-_t & n^{lk} + 1 < L^M \\
 *     ( \mathrm{dn} , n^{lk} , k ) & [ -\Delta^-_t , 0 ] &
 *       ( \mathrm{st} , \tau^M - 1 , k - 1 ) & \Pi_{k-1} & c^-_t &
 *   \end{array} \tag{1}
 * \f]
 * where the moves that start a modulation need \f$ n^M < N^M \f$. Move (i) is
 * the single move that replaces (ii) and (iii) when the direction does not
 * matter (\f$ L^M = 1 \f$, no down cost, no deep decrease and no band); when
 * the output is banded (ii) needs \f$ k < 3 \f$ and (iii) needs
 * \f$ k > 1 \f$, and (iv) needs \f$ L^M > 1 \f$. In the first row the window
 * is the stability ramp, intersected with the ramp of the step; a window
 * \f$ \{ \Delta^+_t \} \f$ is a step of exactly the ramp. When \f$ t + 1 \f$
 * begins a new day, the counters of the landing label are reset to 0. When
 * the deep decreases exist, a move whose window reaches a decrease of
 * \f$ \tilde\Delta_t \f$ is split into three moves. The first one is
 * restricted to a decrease of at least \f$ \tilde\Delta_t \f$ and to
 * \f$ p^{ac}_t \leq \tilde p_t \f$; it is a deep decrease, and therefore it
 * costs \f$ c^{dd}_t \f$ more, needs \f$ n^{dd} < N^{dd} \f$ and adds 1 to
 * \f$ n^{dd} \f$. The second one is restricted to a decrease of at least
 * \f$ \tilde\Delta_t \f$ and to \f$ p^{ac}_t \geq \tilde p_t \f$, and the
 * third one to a decrease of at most \f$ \tilde\Delta_t \f$. At their common
 * boundary the cheaper move, which is not a deep decrease since
 * \f$ c^{dd}_t \geq 0 \f$, prevails, as the strict inequalities of the deep
 * decrease require. A shut-down is allowed only from mode
 * \f$ \mathrm{st} \f$, and it drops the band. An idle instant decreases
 * \f$ n^{lk} \f$ by one, and the counters are reset when a new day begins.
 * Finally, a start-up at \f$ t \f$ from the off-state
 * \f$ ( n^{lk} , n^M , n^{dd} , n^{su} ) \f$ needs \f$ n^{su} < N^{su} \f$
 * and lands in \f$ ( \mathrm{st} , \max\{ n^{lk-} , \tau^v - 1 \} , k ) \f$
 * with \f$ n^{su} + 1 \f$, one state for each band \f$ k \f$ with the range
 * \f$ \Pi_k \f$ of the start-up power. Initially, the label is \f$ (
 * \mathrm{st} , \max\{ \tau^M - \tau^M_0 , 0 \} , k_0 , 0 , 0 , 0 ) \f$, with
 * \f$ k_0 \f$ the lowest band that contains InitialPower. Note that the
 * lockout runs in calendar time, i.e., also along the off instants, and
 * therefore it encodes the window of \f$ \tau^M \f$ instants between two
 * modulations, the initial condition and the stability after a start-up. This
 * is the stability of NuclearUnitBlock, i.e., the rows (7) if \f$ L^M = 1 \f$
 * and (15) otherwise, which are written on \f$ m_t \f$ only, and hence across
 * the off instants as well. For \f$ L^M > 1 \f$ and \f$ \tau^M = 2 \f$ there
 * is no row (15), and none is needed, since, by definition of the end of a
 * modulation, a step right after the last one would continue it. Any label is
 * accepted at \f$ T - 1 \f$. Hence, a modulation that the horizon cuts is cut
 * at \f$ T - 1 \f$ after a step of exactly the ramp, still in the band it
 * started from, with at most \f$ L^M - 1 \f$ steps; as in the rows of
 * NuclearUnitBlock, it cannot have started upwards from the highest band nor
 * downwards from the lowest one (and, as there, it is not required that it
 * could land within \f$ L^M \f$ steps after the horizon). A start-up resets
 * the counters of the day of its landing label when \f$ t + 1 \f$ begins a
 * new day, as the moves do.
 *
 * <b>Domination.</b> A label is at least as good as another one [see
 * label_dominates()] if it has the same mode and the same band, and a lockout
 * (or number of steps) and counters not larger. By (1), such a label allows
 * any sequence of moves, start-ups and shut-downs of the other one, landing
 * again in labels at least as good, as the pruning of ThermalUnitExtDPSolver
 * requires. This pruning is applied with the trimming of the domains [see
 * trim_domination()].
 *
 * <b>Exactness.</b> Move by move, (1) and the rules on the start-ups, the
 * shut-downs and the idle instants are the rows of NuclearUnitBlock at a
 * fixed schedule, and the value functions of the states are exact convex
 * piecewise quadratic functions. Hence, the value of the dynamic program is
 * the optimal value of the rows of NuclearUnitBlock, with the exceptions of
 * ThermalUnitDPSolverBase.h (a ReferenceSchedule, which it refuses, and the
 * parts of the transitions computed by interpolation; the shut-down at 0 is
 * not one, since a NuclearUnitBlock has the ramps). Furthermore, the rule on
 * the extreme bands at the end of the horizon is also written as rows of
 * NuclearUnitBlock, and therefore its abstract representation and this solver
 * have the same optimal value there too. As for the fixings, the fixed
 * Variable of the commitment and of the design are honored as in
 * ThermalUnitExtDPSolver, while those of \f$ m_t \f$ and \f$ d_t \f$, and
 * \f$ \delta_t \f$ fixed to 0, are honored by generating only the moves that
 * agree with them. A \f$ \delta_t \f$ fixed to 1, and a fixed
 * \f$ \delta'_t \f$, \f$ \delta''_t \f$, \f$ s^M_t \f$, \f$ e_t \f$ or
 * \f$ b^k_t \f$ (save the fixings to 0 at \f$ t = 0 \f$ of a unit off before
 * the horizon), are instead refused [see load_fixings()]: the labels do
 * not read the latter, and in the rows the former is the cost and the
 * count of a deep
 * decrease at \f$ t \f$ whatever the unit does there, which is not a move of
 * this program. The spinning reserves, the reactive power and the design
 * variable are handled as for a thermal unit; in particular, the ramp terms
 * of the reserve band are those of the ramps of the step also at a stable
 * instant, since a change of the output due to the reserve is not a
 * modulation. Finally, the schedule is recovered with \f$ m_t = 1 \f$ exactly
 * at the on instants reached by a move with modulation and \f$ d_t = 1 \f$ at
 * those reached by a downward one.
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

#ifndef __NuclearUnitExtDPSolver
 #define __NuclearUnitExtDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitExtDPSolver.h"

#include "NuclearUnitBlock.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS NuclearUnitExtDPSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/// dynamic programming Solver of a NuclearUnitBlock
/** The NuclearUnitExtDPSolver is the ThermalUnitExtDPSolver whose states are
 * labeled by the mode, the lockout, the band and the counters of the day of
 * the unit (cf. the file comment for the model and the algorithm). */

class NuclearUnitExtDPSolver : public ThermalUnitExtDPSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 NuclearUnitExtDPSolver( void ) : ThermalUnitExtDPSolver() {}

 ~NuclearUnitExtDPSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// sets the (NuclearUnitBlock) that the Solver has to solve
 void set_Block( Block * block ) override;

 /// writes the current solution, the modulation included, into the Block
 void get_var_solution( Configuration * solc ) override;

/*--------------------------------------------------------------------------*/
 /// packs the schedule, the modulation included, into a Solution
 /** Returns a NuclearUnitBlockSolution, i.e., what
  * ThermalUnitExtDPSolver::get_Solution() returns plus the modulation and
  * downward indicators of the optimal schedule. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// read the base parameters and, on top, the modulation data
 void load_parameters( void ) override;

 /// intercept the nuclear-specific Modifications (else delegate to base)
 bool guts_of_process_modifications( const p_Mod mod ) override;

 /// read the fixed Variable of the unit, those of the rules comprised
 /** On top of what ThermalUnitExtDPSolver::load_fixings() does, reads the
  * fixed Variable of the operating rules, i.e., the modulation, its
  * direction and the deep decrease: where one of them is fixed, only the
  * moves that agree with its value are generated [see on_moves()], so that
  * the DP solves the problem with the fixings as the rows of
  * NuclearUnitBlock do. The fixings to 0 with which the NuclearUnitBlock
  * encodes the initial conditions (the instants \f$ t < \tau^M - \tau^M_0
  * \f$, and those before the first free commitment of a unit initially off)
  * are read like any other, the label forbidding those moves anyway. A
  * modulation or a direction fixed to 1 at \f$ t \f$ forces the unit on at
  * \f$ t \f$ and forbids a start-up there. A deep decrease fixed to 1 is
  * refused with std::logic_error, since in the rows it is a cost and a
  * count at \f$ t \f$ whatever the unit does there (no row bounds
  * \f$ \delta_t \f$ from above, save (32) of the tight rules), while the
  * DP would keep only the moves that are deep decreases, a stricter
  * problem. A fixed value that is neither 0 nor 1 throws. A fixed Variable
  * that the labels do not read, i.e., \f$ \delta'_t \f$,
  * \f$ \delta''_t \f$, \f$ s^M_t \f$, \f$ e_t \f$ or \f$ b^k_t \f$, is
  * refused with std::logic_error as well, save the fixings to 0 at
  * \f$ t = 0 \f$ of a unit off before the horizon. */
 void load_fixings( void ) override;

 /// the groups of the base solver plus those of the operating rules
 bool reads_group( const std::string & name ) const override;

/*--------------------------------------------------------------------------*/
 // the labels of the states [see the file comment]

 Index on_labels( void ) const override {
  return( f_ncore * f_ncount );
  }

 Index off_labels( void ) const override {
  return( ( lockout_max() + 1 ) * f_ncount );
  }

 Index init_label( void ) const override;

 void on_moves( Index t , Index lab ,
                std::vector< OnMove > & mv ) const override;

 Index shut_label( Index t , Index lab ) const override;

 Index idle_label( Index t , Index e , Index k ) const override;

 Index start_label( Index t , Index e ) const override;

 void start_labels( Index t , Index e ,
                    std::vector< std::pair< Index ,
                                 std::pair< double , double > > > & ls )
  const override;

 bool label_dominates( Index a , Index b ) const override;

 bool trim_domination( void ) const override { return( true ); }

/*--------------------------------------------------------------------------*/
 /// the lockout a modulation leaves behind, \f$ \tau^M - 1 \f$
 Index mod_lockout( void ) const {
  return( std::max( f_mod_interval , Index( 2 ) ) - 1 );
  }

/*--------------------------------------------------------------------------*/
 /// the largest lockout a label may carry,
 /// \f$ \max\{ \tau^M , \tau^v \} - 1 \f$
 /** The largest value the lockout of a label can take: \f$ \tau^M - 1 \f$
  * after the end of a modulation [see mod_lockout()] and \f$ \tau^v - 1 \f$
  * after a start-up, \f$ \tau^v \f$ being the stability that follows one
  * [see NuclearUnitBlock::get_stability_after_start_up()]. It is the range
  * of the lockout in the encoding of the labels, and nothing else. */

 Index lockout_max( void ) const {
  return( std::max( mod_lockout() + 1 , f_stab_start ) - 1 );
  }

 /// the band the output p belongs to, the lowest one at a breakpoint, 0 if
 /// the output is not banded (the bands are numbered 0, 1, 2 here and
 /// 1, 2, 3 in the file comment)
 Index band_of( double p ) const {
  if( f_bands.empty() )
   return( 0 );
  return( ( p <= f_bands[ 0 ] ) ? 0 : ( ( p <= f_bands[ 1 ] ) ? 1 : 2 ) );
  }

/*--------------------------------------------------------------------------*/
 /// the day of the time instant t
 Index day( Index t ) const { return( f_day_length ? t / f_day_length : 0 ); }

 /// a label, decoded: the mode \f$ \varpi \f$ (0 stable, 1 up, 2 down),
 /// the lockout or the steps \f$ n^{lk} \f$, the counters of the day
 /// \f$ n^M \f$, \f$ n^{dd} \f$, \f$ n^{su} \f$ and the band \f$ k \f$
 struct Label {
  int mode;
  Index lk;
  Index c , a , s;
  Index b{};   ///< the band of the output, or the one a modulation left
  };

 /// encode an on-label
 Index on_code( const Label & l ) const {
  const Index core = ( l.mode == 0 ) ? l.lk :
   lockout_max() + l.lk + ( l.mode == 2 ? f_max_mod_length - 1 : 0 );
  return( core + f_ncore * count_code( l ) );
  }

 /// encode an off-label (always stable)
 Index off_code( const Label & l ) const {
  return( l.lk + ( lockout_max() + 1 ) * count_code( l ) );
  }

 /// decode an on-label
 Label on_label( Index lab ) const;

 /// decode an off-label
 Label off_label( Index e ) const;

 /// the counters part of a label
 Index count_code( const Label & l ) const {
  return( l.b + f_nband * ( l.c + f_nc * ( l.a + f_na * l.s ) ) );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 Index f_mod_interval{ 2 };            ///< \f$ \tau^M \f$ (ModulationTime)
 int f_init_modulation{ 2 };           ///< \f$ \tau^M_0 \f$ (InitModulation)
 std::vector< double > mod_ramp_up;    ///< \f$ \Delta^{M+}_t \f$
 std::vector< double > mod_ramp_down;  ///< \f$ \Delta^{M-}_t \f$

 Index f_max_mod_length{ 1 };          ///< \f$ L^M \f$
 int f_mod_per_day{ -1 };              ///< \f$ N^M \f$, -1 if unlimited
 int f_deep_per_day{ -1 };             ///< \f$ N^{dd} \f$, -1 if unlimited
 int f_starts_per_day{ -1 };           ///< \f$ N^{su} \f$, -1 if unlimited
 Index f_day_length{ 0 };              ///< \f$ T^{day} \f$, 0 = horizon
 bool f_direction{ false };            ///< true if the direction matters
 std::vector< double > f_down_cost;    ///< \f$ c^-_t \f$ (empty = 0)
 std::vector< double > f_deep_thr;     ///< \f$ \tilde p_t \f$ (empty = none)
 std::vector< double > f_deep_grad;    ///< \f$ \tilde\Delta_t \f$
 std::vector< double > f_deep_cost;    ///< \f$ c^{dd}_t \f$ (empty = 0)

 /// \f$ \tau^v \f$, the instants of stability that begin with a start-up
 Index f_stab_start{};

 /// the two breakpoints that split the output into bands, empty if there
 /// are no bands [see NuclearUnitBlock::get_power_bands()]
 std::vector< double > f_bands;

 /// the number of bands: 3 with the breakpoints, 1 without
 Index f_nband{ 1 };

 /// the fixed Variable \f$ m_t \f$, \f$ d_t \f$ and \f$ \delta_t \f$ of the
 /// operating rules [see load_fixings()]: -1 where free, 0 or 1 where
 /// fixed, empty if none of them is fixed
 std::vector< signed char > f_fix_mod;
 std::vector< signed char > f_fix_down;  ///< as f_fix_mod, for \f$ d_t \f$
 std::vector< signed char > f_fix_deep;  ///< the same, for \f$ \delta_t \f$

 /// true if a Variable of the operating rules is fixed to 1 at t, which
 /// needs the unit on at t and not starting up at t
 bool on_forced( Index t ) const {
  auto one = [ t ]( const std::vector< signed char > & f ) {
   return( ( ! f.empty() ) && ( f[ t ] == 1 ) );
   };
  return( one( f_fix_mod ) || one( f_fix_down ) || one( f_fix_deep ) );
  }

 Index f_ncore{ 2 };   ///< number of (mode, lockout or steps) pairs
 Index f_nc{ 1 };      ///< range of the counter of the modulations
 Index f_na{ 1 };      ///< range of the counter of the deep decreases
 Index f_nv{ 1 };      ///< range of the counter of the start-ups
 Index f_ncount{ 1 };  ///< f_nband * f_nc * f_na * f_nv

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( NuclearUnitExtDPSolver ) )

/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* NuclearUnitExtDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitExtDPSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
