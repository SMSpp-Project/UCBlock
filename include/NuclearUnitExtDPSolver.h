/*--------------------------------------------------------------------------*/
/*--------------------- File NuclearUnitExtDPSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the NuclearUnitExtDPSolver class, a Solver for the
 * NuclearUnitBlock that solves the single-Unit Commitment (1UC) problem of a
 * nuclear unit, i.e., of a thermal unit subject to modulation constraints,
 * by Dynamic Programming.
 *
 * <b>The model.</b> On top of the whole ThermalUnitBlock model, a
 * NuclearUnitBlock has a binary modulation indicator \f$ m_t \f$ per time
 * instant and the constraints
 * \f[
 * \begin{array}{ll}
 *  p_t - p_{t-1} \leq \Delta^{M+}_t u_{t-1} + ( \Delta^+_t - \Delta^{M+}_t )
 *  m_t + SU_t v_t & \mbox{(modulation ramp-up)} \\
 *  p_{t-1} - p_t \leq \Delta^{M-}_t u_t + ( \Delta^-_t - \Delta^{M-}_t ) m_t
 *  + SD_t w_t & \mbox{(modulation ramp-down)} \\
 *  m_t \leq u_t \, , \quad m_t + v_t \leq 1 & \\
 *  \sum_{h = \max\{ 0 , t - \tau^M + 1 \}}^t m_h \leq 1 &
 *  \mbox{(modulation window)}
 * \end{array}
 * \f]
 * with \f$ v_t \f$ and \f$ w_t \f$ the start-up and shut-down indicators,
 * \f$ 0 \leq \Delta^{M\pm}_t \leq \Delta^\pm_t \f$ the modulation ramps
 * (ModulationDeltaRampUp / ModulationDeltaRampDown), \f$ \tau^M \geq 2 \f$
 * the modulation interval (ModulationTime) and the initial condition that
 * fixes \f$ m_t = 0 \f$ for \f$ 0 \leq t < \tau^M - \kappa^0 \f$, where the
 * unit last modulated \f$ \kappa^0 \geq 1 \f$ instants before 0
 * (InitModulation). That is, along a continuous on-run the power moves by
 * at most the modulation ramp \f$ \Delta^{M\pm}_t \f$, unless a modulation
 * is performed, which allows the full ramp \f$ \Delta^\pm_t \f$; two
 * modulations are at least \f$ \tau^M \f$ instants apart; start-ups and
 * shut-downs are not modulations and keep their caps \f$ SU_t \f$ and
 * \f$ SD_t \f$. Since the ramp constraints of the ThermalUnitBlock are still
 * there, the window of the move from \f$ t - 1 \f$ to \f$ t \f$ of a unit
 * that stays on is, depending on \f$ m_t \f$,
 * \f[
 *  -\min\{ \Delta^-_{t-1} , \widehat{\Delta}^-_t( m_t ) \} \leq p_t - p_{t-1}
 *  \leq \min\{ \Delta^+_{t-1} , \widehat{\Delta}^+_t( m_t ) \}
 *  \, , \quad \widehat{\Delta}^\pm_t( m ) = m \Delta^\pm_t + ( 1 - m )
 *  \Delta^{M\pm}_t
 * \f]
 * (the ramps of the ThermalUnitBlock being indexed by the first instant of
 * the step and those of the modulation by the second one), while the
 * reserve held at \f$ t \f$ keeps the deliverability of the full ramp
 * \f$ \Delta^\pm_{t-1} \f$: a change of the output due to the reserve is
 * not a modulation.
 *
 * <b>The lockout counter.</b> The window constraint couples the \f$ m_t \f$
 * across time, so (unlike the reserves) it cannot be folded into a
 * per-period cost; since it says that any two modulations are at least
 * \f$ \tau^M \f$ instants apart, it is enough to remember how long ago the
 * unit last modulated, capped at \f$ \tau^M - 1 \f$. The *lockout*
 * \f$ \ell \in \{ 0 , \ldots , \tau^M - 1 \} \f$ entering \f$ t \f$ is the
 * number of instants, from \f$ t \f$ on, in which a modulation is still
 * forbidden, and its dynamics is
 * \f[
 *  \ell_{t+1} = \left\{ \begin{array}{ll}
 *   \tau^M - 1 & \mbox{if } m_t = 1 \mbox{ (allowed only if } \ell_t = 0
 *   \mbox{)} \\
 *   \max\{ \ell_t - 1 , 0 \} & \mbox{if } m_t = 0
 *  \end{array} \right.
 * \f]
 * in calendar time, i.e., also along the off instants (where \f$ m_t = 0
 * \f$), starting from \f$ \ell_0 = \max\{ \tau^M - \kappa^0 , 0 \} \f$. The
 * lockout is the label of the states of the DP of ThermalUnitExtDPSolver
 * [see the label methods there], which carries it through the value
 * functions \f$ F^{\tau,\ell}_t \f$ of the on-states and through the
 * ready-to-restart values \f$ c^{rdy}_t( \ell ) \f$ of the off-states. An
 * on-state with lockout \f$ \ell \f$ continues by the move without
 * modulation, of window \f$ \widehat{\Delta}^\pm_t( 0 ) \f$ and landing
 * lockout \f$ \max\{ \ell - 1 , 0 \} \f$, and, if \f$ \ell = 0 \f$, by the
 * move with modulation, of window \f$ \widehat{\Delta}^\pm_t( 1 ) \f$ and
 * landing lockout \f$ \tau^M - 1 \f$; a shut-down carries the lockout into
 * the off-states, the idle instants decrease it and a restart (with
 * \f$ m_t = 0 \f$ forced) decreases it once more. A smaller lockout is at
 * least as good as a larger one, since it allows the same moves or more,
 * which is what the domination pruning uses. Every move is a sliding
 * minimum over its own window, so the value functions stay convex piecewise
 * quadratic, the states with different lockouts being kept apart rather
 * than merged into a (non-convex) pointwise minimum, and \f$ m_t = 1 \f$
 * exactly at the instants reached by the move with modulation along the
 * optimal schedule. The state grows by a factor \f$ \tau^M \f$, hence the
 * worst case is \f$ O( n^3 \tau^M ) \f$, and the spinning reserves, the
 * reactive power, the design variable and the fixing of the commitment are
 * handled exactly as for a thermal unit.
 *
 * <b>The operating rules.</b> The above is the case of single-instant
 * modulations and no daily limit; in general [see NuclearUnitBlock] a
 * modulation lasts up to \f$ L^M \f$ instants, all but the last at the full
 * ramp, at most \f$ C \f$ modulations, \f$ A \f$ deep decreases and
 * \f$ V \f$ start-ups happen in each day, and the downward steps and the
 * deep decreases have a cost. The label of an on-state entering \f$ t + 1
 * \f$ is then the tuple \f$ ( \mu , \ell , k , c , a , s ) \f$: the mode
 * \f$ \mu \in \{ S , U , D \} \f$ (stable, or in the middle of an upward
 * or downward modulation), the lockout \f$ \ell \f$ (if stable), the number
 * \f$ k \in \{ 1 , \ldots , L^M - 1 \} \f$ of steps already taken (if
 * not), and the numbers \f$ c , a , s \f$ of modulations, deep decreases
 * and start-ups in the day of \f$ t + 1 \f$, reset when it starts a new
 * day; the label of an off-state is \f$ ( \ell , c , a , s ) \f$. The moves
 * out of an on-state at \f$ t - 1 \f$ are, with \f$ \ell^- = \max\{ \ell -
 * 1 , 0 \} \f$ and \f$ \tau^M = B + 1 \f$:
 * \f[
 * \begin{array}{llll}
 *  \mbox{from} & \mbox{window of } p_t - p_{t-1} & \mbox{to} &
 *  \mbox{when} \\
 *  ( S , \ell ) & [ -\Delta^{M-}_t , \Delta^{M+}_t ] & ( S , \ell^- ) & \\
 *  ( S , 0 ) & [ 0 , \Delta^+_t ] \mbox{ or } [ -\Delta^-_t , 0 ] &
 *  ( S , B ) , c + 1 & c < C \\
 *  ( S , 0 ) & \{ \Delta^+_t \} \mbox{ or } \{ -\Delta^-_t \} &
 *  ( U , 1 ) \mbox{ or } ( D , 1 ) , c + 1 & c < C \, , \, L^M > 1 \\
 *  ( U , k ) & \{ \Delta^+_t \} & ( U , k + 1 ) & k < L^M - 1 \\
 *  ( U , k ) & [ 0 , \Delta^+_t ] & ( S , B ) & \\
 *  ( D , k ) & \{ -\Delta^-_t \} & ( D , k + 1 ) & k < L^M - 1 \\
 *  ( D , k ) & [ -\Delta^-_t , 0 ] & ( S , B ) &
 * \end{array}
 * \f]
 * (each window intersected with the ramp \f$ [ -\Delta^-_{t-1} ,
 * \Delta^+_{t-1} ] \f$ of the ThermalUnitBlock, a full-ramp step being
 * impossible if the latter is tighter), each downward move costing
 * \f$ c^-_t \f$. A move whose window reaches a decrease of
 * \f$ \tilde{\Delta}_t \f$ is split into the part with a decrease of at
 * least \f$ \tilde{\Delta}_t \f$ and a landing power of at most
 * \f$ \tilde{p}_t \f$ (a deep decrease: cost \f$ c^d_t \f$, \f$ a + 1
 * \leq A \f$), the part with the same decrease and a landing power of at
 * least \f$ \tilde{p}_t \f$, and the part with a decrease of at most
 * \f$ \tilde{\Delta}_t \f$, the three overlapping only at their boundary,
 * where the cheaper (not deep) one prevails. Only a stable unit can shut
 * down, the idle instants decrease the lockout, and a restart needs
 * \f$ s < V \f$, increases \f$ s \f$ and gives \f$ ( S , \ell^- ) \f$.
 * A label dominates another one with the same mode and at least as large
 * lockout, number of steps and counters. With single-instant modulations
 * and no down cost nor deep decrease the two modulation windows are merged
 * into \f$ [ -\Delta^-_t , \Delta^+_t ] \f$, which gives back the lockout
 * counter above.
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
/// DP solver for the 1UC of a nuclear unit (thermal + modulation)
/** The NuclearUnitExtDPSolver is the ThermalUnitExtDPSolver whose states are
 * labelled by the modulation lockout; see the file comment for the model and
 * the algorithm. */

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
  * ThermalUnitExtDPSolver::get_Solution() returns plus the modulation
  * indicators of the optimal schedule. */

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
  * the DP solves the problem with the fixings, exactly as a MILP solver
  * would. The fixings to 0 with which the NuclearUnitBlock encodes the
  * initial conditions (the instants \f$ t < \tau^M - \kappa^0 \f$, and
  * those before the first free commitment of a unit initially off) are
  * read like any other, the label forbidding those moves anyway. A fixed
  * value that is neither 0 nor 1 throws. */
 void load_fixings( void ) override;

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
 /// the largest lockout a label may carry, \f$ \max\{ \tau^M , A \} - 1 \f$
 /** The largest value the lockout of a label can take: \f$ \tau^M - 1 \f$
  * after the end of a modulation [see mod_lockout()] and \f$ A - 1 \f$
  * after a start-up, A being the stability that follows one [see
  * NuclearUnitBlock::get_stability_after_start_up()]. It is the range of
  * the lockout in the encoding of the labels, and nothing else. */

 Index lockout_max( void ) const {
  return( std::max( mod_lockout() + 1 , f_stab_start ) - 1 );
  }

 /// the band the output p belongs to, 0 if the output is not banded
 Index band_of( double p ) const {
  if( f_bands.empty() )
   return( 0 );
  return( ( p <= f_bands[ 0 ] ) ? 0 : ( ( p <= f_bands[ 1 ] ) ? 1 : 2 ) );
  }

/*--------------------------------------------------------------------------*/
 /// the day of the time instant t
 Index day( Index t ) const { return( f_day_length ? t / f_day_length : 0 ); }

 /// a label, decoded: mode (0 stable, 1 up, 2 down), lockout or steps, and
 /// the counters of the day
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

 Index f_mod_interval{ 2 };            ///< modulation interval tau^M
 int f_init_modulation{ 2 };           ///< InitModulation kappa^0
 std::vector< double > mod_ramp_up;    ///< modulation Delta^M+ per instant
 std::vector< double > mod_ramp_down;  ///< modulation Delta^M- per instant

 Index f_max_mod_length{ 1 };          ///< L^M
 int f_mod_per_day{ -1 };              ///< C, -1 if unlimited
 int f_deep_per_day{ -1 };             ///< A, -1 if unlimited
 int f_starts_per_day{ -1 };           ///< V, -1 if unlimited
 Index f_day_length{ 0 };              ///< instants of a day, 0 = horizon
 bool f_direction{ false };            ///< true if the direction matters
 std::vector< double > f_down_cost;    ///< c^-_t (empty = 0)
 std::vector< double > f_deep_thr;     ///< deep threshold (empty = none)
 std::vector< double > f_deep_grad;    ///< deep gradient
 std::vector< double > f_deep_cost;    ///< c^d_t (empty = 0)

 /// the fixed Variable of the operating rules [see load_fixings()]: -1
 /// where free, 0 or 1 where fixed, empty if none of them is fixed
 /// the instants of stability that follow a start-up
 Index f_stab_start{};

 /// the two breakpoints that split the output into bands, empty if there
 /// are no bands [see NuclearUnitBlock::get_power_bands()]
 std::vector< double > f_bands;

 /// the number of bands: 3 with the breakpoints, 1 without
 Index f_nband{ 1 };

 std::vector< signed char > f_fix_mod;
 std::vector< signed char > f_fix_down;
 std::vector< signed char > f_fix_deep;

 Index f_ncore{ 2 };   ///< number of (mode, lockout or steps) pairs
 Index f_nc{ 1 };      ///< range of the counter of the modulations
 Index f_na{ 1 };      ///< range of the counter of the deep decreases
 Index f_nv{ 1 };      ///< range of the counter of the start-ups
 Index f_ncount{ 1 };  ///< f_nc * f_na * f_nv

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
