/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitDPSolverBase.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolverBase class, the common base of the
 * two dynamic programming Solvers of the single-unit commitment problem of a
 * ThermalUnitBlock. These are ThermalUnitDPSolver, which enumerates the
 * on-intervals of the unit on an acyclic graph and prices each of them by an
 * economic dispatch, and ThermalUnitExtDPSolver, which carries one value
 * function per length of the current on-run (and per label, for a derived
 * unit).
 *
 * Both solve the problem given by the rows of ThermalUnitBlock [see
 * ThermalUnitBlock::generate_abstract_constraints()], which we recall in the
 * notation of that class. Over the instants \f$ t \in \mathcal{T} = \{ 0 ,
 * \ldots , T - 1 \} \f$, the commitment \f$ u_t \f$ obeys the minimum up and
 * down times \f$ \tau^+ \f$ and \f$ \tau^- \f$ from the initial state
 * (\f$ \tau_0 \f$, \f$ t_0 \f$, \f$ p_{-1} \f$). Also, the active power
 * \f$ p^{ac}_t \f$ of an on instant lies in \f$ [ P^{mn}_t , P^{mx}_t ] \f$,
 * it is at most \f$ P^{su}_t \f$ at a start-up at \f$ t \f$ and at most
 * \f$ P^{sd}_{t+1} \f$ when the unit is off at \f$ t + 1 \f$, and it moves
 * from \f$ t - 1 \f$ to \f$ t \f$ by at most \f$ \Delta^+_t \f$ upwards and
 * \f$ \Delta^-_t \f$ downwards while the unit stays on (the step from the
 * initial state into \f$ t = 0 \f$ is indexed by 0). The cost of an on
 * instant is
 * \f[
 *   f_t( p ) = a_t p^2 + b_t p + c_t \tag{1}
 * \f]
 * (QuadTerm, LinearTerm and ConstTerm, \f$ a_t \geq 0 \f$, while \f$ b_t \f$
 * may carry a Lagrangian price of any sign), a start-up at \f$ t \f$ costs
 * \f$ c^{su}_t \f$ and a shut-down at \f$ t \f$, i.e., the unit on at
 * \f$ t - 1 \f$ and off at \f$ t \f$, costs \f$ c^{sd}_t \f$. We use the
 * operational bounds of the rows of ThermalUnitBlock,
 * \f$ P^{mx}_t = \chi_t \hat P^{mx}_t \f$ and
 * \f$ P^{mn}_t = \hat P^{mn}_t \f$ (0 if \f$ \chi_t = 0 \f$), with the
 * availability \f$ \chi_t \f$ and the nominal MaxPower \f$ \hat P^{mx}_t \f$
 * and MinPower \f$ \hat P^{mn}_t \f$. When DeltaRampUp (DeltaRampDown) is not
 * given, ThermalUnitBlock has no ramp row of that kind, and
 * \f$ \Delta^+_t \f$ (\f$ \Delta^-_t \f$) is taken to be the largest of
 * InitialPower \f$ p_{-1} \f$ and of the \f$ \hat P^{mx}_t \f$ over the
 * horizon, which no move exceeds, also the one from \f$ p_{-1} \f$ to the
 * power at 0; hence, the moves and the reserves are free, as in the rows,
 * and a change of InitialPower makes the solvers compute it again. If
 * "FixToMaximum" is positive, the rows
 * \f$ p^{ac}_t \geq P^{mx}_t \f$ make the output of an on instant
 * \f$ P^{mx}_t \f$ and the unit on wherever \f$ P^{mx}_t > 0 \f$. In this
 * case the solvers take \f$ P^{mn}_t = P^{mx}_t \f$ and fix the unit on at
 * those instants [see force_on_fixed_to_maximum()]; this leaves the reserves
 * unchanged, since their room above the output is 0 in both cases.
 *
 * This class holds what the two solvers share, i.e., the data loaded from the
 * ThermalUnitBlock [see load_common_parameters()], the convex piecewise
 * quadratic functions of one variable in which the value functions are stored
 * (PieceQuad, PQFun) with the operations on them, and the model of the
 * spinning reserves. We now describe the latter. At an on instant \f$ t \f$
 * with active power \f$ p \f$ the primary and secondary reserves
 * \f$ p^{pr}_t \f$ and \f$ p^{sc}_t \f$ enter only the objective of the unit,
 * with the costs \f$ c^{pr}_t \f$ and \f$ c^{sc}_t \f$. A reserve exists only
 * if the unit has its Variable, i.e., if the enclosing UCBlock requires it
 * and the unit has the fraction [see
 * ThermalUnitBlock::has_primary_reserve()]; otherwise, its fraction and cost
 * are ignored, as the rows of ThermalUnitBlock have neither. Both reserves
 * are bounded by the fractions \f$ \rho^{pr}_t \f$ and \f$ \rho^{sc}_t \f$ of
 * \f$ p \f$ and, together, by a band \f$ H \f$; hence, their contribution to
 * the cost of the instant is the value of the linear program
 * \f[
 *   g_t( p , H ) = \min \{ \, c^{pr}_t r^{pr} + c^{sc}_t r^{sc} \, : \,
 *     0 \leq r^{pr} \leq \rho^{pr}_t p \, , \;
 *     0 \leq r^{sc} \leq \rho^{sc}_t p \, , \;
 *     r^{pr} + r^{sc} \leq H \, \} \tag{2}
 * \f]
 * (0 if \f$ H \leq 0 \f$, a value that the solvers never use; cf. (4) below).
 * This linear program is solved by filling the band with the reserve of the
 * most negative cost first, and then with the other one if its cost is
 * negative too. In particular, \f$ g_t \equiv 0 \f$ (and the reserves are 0)
 * when neither cost is negative, as it happens when the unit is solved alone
 * with nonnegative reserve costs. The band, in turn, is
 * \f[
 *   \mathcal{H}_t( p , q ) = \min \{ \, p - P^{mn}_t \, , \; K_t - p \, , \;
 *     \Delta^+_t - ( p - q ) \, , \; \Delta^-_t + ( p - q ) \, \}
 *   \tag{3}
 * \f]
 * where \f$ q \f$ is the active power at \f$ t - 1 \f$, and the last two
 * terms are present only if the unit is on at \f$ t - 1 \f$ and
 * \f$ t \geq 1 \f$. The cap \f$ K_t \f$ is \f$ P^{mx}_t \f$ at an instant
 * inside an on-run, \f$ P^{su}_t \f$ at the first instant of an on-run
 * started at \f$ t \f$, and the smaller of that value and
 * \f$ P^{sd}_{t+1} \f$ when the unit is off at \f$ t + 1 \f$. These are
 * exactly the rows of the reserves of ThermalUnitBlock (the room below the
 * output, the room above it with the start-up and shut-down limits, the
 * deliverability within the ramp left by the scheduled move, and the
 * fractions), solved for the reserves at a fixed schedule. At \f$ t = 0 \f$
 * with the unit on before the horizon, \f$ q \f$ is InitialPower
 * \f$ p_{-1} \f$, as in the deliverability rows of ThermalUnitBlock at 0;
 * hence, the two ramp terms of (3) are a floor and a cap of the power, and
 * the band is \f$ \min\{ p - \max\{ P^{mn}_0 , p_{-1} - \Delta^-_0 \} ,
 * \min\{ K_0 , p_{-1} + \Delta^+_0 \} - p \} \f$ [see
 * initial_reserve_discount()]. Finally, the costs \f$ c^{pr}_t \f$ and
 * \f$ c^{sc}_t \f$ are PrimarySpinningReserveCost and
 * SecondarySpinningReserveCost (0 when a vector is not given), as in
 * ThermalUnitBlock::generate_objective(). Its Configuration only determines
 * whether a reserve with cost 0 has a term in the Objective, which makes no
 * difference for the solvers.
 *
 * The transition of an on-run from \f$ t - 1 \f$ to \f$ t \f$ is
 * \f[
 *   G( p ) = \min \{ \, F( q ) + g_t( p , \mathcal{H}_t( p , q ) ) \, : \,
 *     p - \omega^+ \leq q \leq p + \omega^- \, \} \tag{4}
 * \f]
 * for a convex value function \f$ F \f$ of the power at \f$ t - 1 \f$ and a
 * window \f$ [ -\omega^- , \omega^+ ] \f$ of the move \f$ p - q \f$, with
 * \f$ \omega^\pm \leq \Delta^\pm_t \f$ and \f$ p \in [ P^{mn}_t , K_t ] \f$.
 * On this domain every term of (3) is nonnegative, and therefore
 * \f$ \mathcal{H}_t \geq 0 \f$ and the convention for \f$ H \leq 0 \f$ in (2)
 * never acts. Thus, \f$ ( q , p ) \mapsto g_t( p , \mathcal{H}_t( p , q )
 * ) \f$ is jointly convex there, since it is the value of the linear program
 * (2) whose right-hand side is affine in \f$ ( p , q ) \f$ (which the
 * convention would break across \f$ H = 0 \f$ when a cost is negative);
 * \f$ G \f$ is then convex and piecewise quadratic. This function is computed
 * exactly by sliding_min() when no reserve cost is negative, and otherwise by
 * sliding_min_corr(). In this case the computation is exact except on the
 * parts of the domain that it covers by linear interpolation of exact values
 * of \f$ G \f$; these are secants of the convex \f$ G \f$, i.e., an
 * approximation from above within a tolerance of
 * \f$ 10^{-7} \max\{ 1 , | G | \} \f$ (or on intervals shorter than
 * \f$ 10^{-5} \f$) [see sliding_min_corr()]. For a
 * thermal unit the window is the ramp \f$ [ -\Delta^-_t , \Delta^+_t ] \f$ of
 * the step, while for a derived unit it may be a narrower one (a
 * NuclearUnitBlock that does not modulate). However, the last two terms of
 * (3) keep the ramp of the step: a change of the output due to the reserve is
 * not a modulation.
 *
 * Furthermore, the values reported by the solvers are multiplied by the scale
 * factor \f$ \sigma \f$ of the unit [see scaled_value()]. The solvers find
 * the optimum of the model above, i.e., of the rows of ThermalUnitBlock, with
 * the following exceptions:
 *
 * - a unit with a ReferenceSchedule, whose Objective has a term these
 *   solvers do not represent, is refused by load_common_parameters();
 *
 * - a fixed Variable other than the commitment and the design (and those
 *   that ThermalUnitBlock fixes for the initial state) is refused by the
 *   load_fixings() of each solver;
 *
 * - the parts of (4) computed by interpolation, above.
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

#ifndef __ThermalUnitDPSolverBase
 #define __ThermalUnitDPSolverBase
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "ThermalUnitBlock.h"

#include <cmath>
#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ThermalUnitDPSolverBase -----------------------*/
/*--------------------------------------------------------------------------*/
/// common base of the two dynamic programming Solvers of a ThermalUnitBlock
/** Holds the data loaded from the ThermalUnitBlock, the convex piecewise
 * quadratic value functions with their operations and the model of the
 * spinning reserves, i.e., (2)-(4) of the file comment. Although it derives
 * from Solver, it is not a complete Solver: the structure of the dynamic
 * program and the interface of the Solver (compute(), get_var_solution(), and
 * so on) are those of ThermalUnitDPSolver and ThermalUnitExtDPSolver. */

class ThermalUnitDPSolverBase : public Solver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

 static constexpr auto TUEDPINF = Inf< double >();  ///< the INF value

 using Index = Block::Index;

 ThermalUnitDPSolverBase( void ) {}

 ~ThermalUnitDPSolverBase() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
 /// the value of \p v for all the copies of the unit
 /** The DP works on the data of one copy of the unit, while its Objective
  * is the scale factor \f$ \sigma \f$ times the cost of one copy [see
  * UnitBlock::scale()]: the values it reports are \f$ \sigma v \f$. A unit
  * with \f$ \sigma = 0 \f$ costs nothing, unless it is infeasible. */

 OFValue scaled_value( OFValue v ) const {
  const auto scale = f_Block ?
   static_cast< const UnitBlock * >( f_Block )->get_scale() : 1.0;
  if( scale == 1 )
   return( v );
  if( scale == 0 )
   return( std::isinf( v ) ? v : 0 );
  return( scale * v );
  }

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// one quadratic piece \f$ \alpha p^2 + \beta p + \gamma \f$ on an
 /// interval of the power \f$ p \f$
 /** The piece of a piecewise quadratic value function: the quadratic
  * \f$ \alpha p^2 + \beta p + \gamma \f$ on the interval
  * \f$ [ \mbox{left} , \mbox{right} ] \f$ of the power. Two consecutive
  * pieces of a PQFun share their common endpoint. */

 struct PieceQuad {
  double alfa;   ///< coefficient of \f$ p^2 \f$ (>= 0 for a convex piece)
  double beta;   ///< coefficient of \f$ p \f$
  double gamma;  ///< additive constant (absorbs accumulated path cost)
  double left;   ///< left endpoint of the piece (inclusive)
  double right;  ///< right endpoint of the piece (inclusive)
 };

 /// a convex piecewise quadratic function on the power axis
 /** Stored as a vector of PieceQuad sorted by their left endpoint, which
  * cover an interval of the power with no gap and no overlap but at the
  * shared endpoints; on it the function is
  *
  * - defined piece by piece, the right endpoint of each piece being the left
  *   endpoint of the next one;
  *
  * - continuous, the two pieces that share an endpoint having the same value
  *   there;
  *
  * - convex, each piece having \f$ \alpha \geq 0 \f$, with a derivative
  *   that is nondecreasing across the endpoints.
  *
  * The function is \f$ +\infty \f$ outside that interval, and the empty
  * vector is the function \f$ +\infty \f$ everywhere, i.e., a state of
  * the dynamic program that is infeasible or not reached. All the
  * operations below keep these properties. */

 using PQFun = std::vector< PieceQuad >;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// load the data shared by both solvers from the ThermalUnitBlock
 /** Reads (under a read-lock) the scalar parameters, the power/ramp/cost
  * bounds, the reactive box, the spinning-reserve participation factors and
  * prices and the design cost into the protected fields below. Each
  * subclass' load_parameters() calls this and then initialises its own
  * DP-specific state. Requires f_Block to be a ThermalUnitBlock (checked in
  * set_Block). */
 void load_common_parameters( void );

 /// adds the instants at which "FixToMaximum" keeps the unit on
 /** If the unit is fixed to its maximum power [see
  * ThermalUnitBlock::is_fixed_to_maximum()], marks as fixed on, in the
  * tables of the fixings, every instant whose maximum power is positive;
  * each load_fixings() calls it after having read the fixed Variable. */
 void force_on_fixed_to_maximum( void );

 /// the name of a group of Variable of an extended formulation with a fixing
 /** Returns the name of the first group of static Variable of the
  * ThermalUnitBlock that the solver does not read [see reads_group()] and
  * that has a fixed element, or the empty string if there is none;
  * ThermalUnitBlock fixes none of them, and the solvers cannot honour such
  * a fixing, which each load_fixings() therefore refuses. */
 std::string fixed_extended_variable( void ) const;

 /// whether load_fixings() reads the group of static Variable \p name
 /** True for the groups of the commitment, the design, the active and
  * reactive power, the reserves, the start-up and the shut-down; a solver
  * of a derived unit adds the groups whose fixings it reads itself. */
 virtual bool reads_group( const std::string & name ) const;

 /// expand a per-instant vector: an empty one is 0, a single value is
 /// repeated over the horizon
 void retrieve_term( std::vector< double > & out ,
                     const std::vector< double > & in ) const;

/*--------------------------------------------------------------------------*/
/*---------------------------- RESERVE MODEL -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// optimal reserves at instant t with active power p and cap \p cap
 /** Solves (2) at \p t with \f$ H = \min\{ p - \mbox{floor} ,
  * \mbox{cap} - p \} \f$, where \p floor is taken as \f$ P^{mn}_t \f$
  * when it is NaN, i.e., the band (3) without the ramp terms, writes the
  * optimal reserves in \p pr and \p sr and returns \f$ g_t( p , H ) \f$,
  * which is 0 unless some reserve cost is negative. */
 double reserve_alloc( Index t , double p , double & pr , double & sr ,
                       double cap ,
                       double floor = std::nan( "" ) ) const;

 /// packs a schedule into a ThermalUnitBlockSolution
 /** Packs the schedule the dynamic programming has found, i.e., the active
  * power \p p, the commitment \p u, the primary and secondary spinning
  * reserve \p pr and \p sr, the reactive power \p q and the value
  * \p design of the dimensioning variable, into a new
  * ThermalUnitBlockSolution [see ThermalUnitBlock.h]. Each of the parts is
  * only saved if the corresponding vector is nonempty, while the active
  * power and the commitment are always there; the Solution is filled
  * directly, so that the ThermalUnitBlock is not written into and no
  * Variable is required to exist. Ownership of the returned object is the
  * caller's. */

 Solution * pack_Solution( const std::vector< double > & p ,
                           const std::vector< double > & u ,
                           const std::vector< double > & pr ,
                           const std::vector< double > & sr ,
                           const std::vector< double > & q ,
                           double design ) const;

/*--------------------------------------------------------------------------*/
 /// true if the ThermalUnitBlock has reactive power data

 bool has_reactive_power( void ) const;

/*--------------------------------------------------------------------------*/
 /// optimal reserves at instant t with active power p and band \p H
 /** Solves (2) at \p t for the given band \p H, writes the optimal
  * reserves in \p pr and \p sr and returns \f$ g_t( p , H ) \f$. */
 double reserve_alloc_band( Index t , double p , double H ,
                            double & pr , double & sr ) const;

 /// the reserve term of an instant with no ramp coupling, as a PQFun
 /** Returns the convex piecewise linear function
  * \f$ p \mapsto g_t( p , \min\{ p - P^{mn}_t , \mbox{cap} - p \} ) \f$
  * on \f$ [ P^{mn}_t , \mbox{cap} ] \f$, i.e., (2) with the band (3)
  * without the ramp terms and with \f$ K_t \f$ equal to \p cap, which is
  * what an instant with no on predecessor contributes (the first instant of
  * an on-run, and instant 0 of a unit on before the horizon). It is
  * nonpositive, and the empty PQFun, i.e., \f$ g_t \equiv 0 \f$, when no
  * reserve cost is negative. The cap \p cap is \f$ P^{mx}_t \f$,
  * \f$ P^{su}_t \f$ or the smaller of those and \f$ P^{sd}_{t+1} \f$ as
  * in (3); \p floor, if not NaN, takes the place of \f$ P^{mn}_t \f$. */
 PQFun build_reserve_discount( Index t , double cap ,
                               double floor = std::nan( "" ) ) const;

 /// the reserve term of instant 0 of a unit on before the horizon
 /** Returns the function \f$ p \mapsto g_0( p , \mathcal{H}_0( p , p_{-1}
  * ) ) \f$ of (2) and (3) on its domain, where \f$ p_{-1} \f$ is InitialPower
  * and the cap \f$ K_0 \f$ is \p cap: with a fixed \f$ q = p_{-1} \f$ the
  * two ramp terms of (3) are a floor and a cap of the power, and the band is
  * \f$ \min\{ p - \max\{ P^{mn}_0 , p_{-1} - \Delta^-_0 \} ,
  * \min\{ K_0 , p_{-1} + \Delta^+_0 \} - p \}
  * \f$ [see build_reserve_discount()]. */
 PQFun initial_reserve_discount( double cap ) const;

 /// the value \f$ g_t( p , H ) \f$ of (2), without the reserves
 /** Returns \f$ g_t( p , H ) \f$, which is 0 when \f$ H \leq 0 \f$ or
  * when no reserve cost is negative; used by sliding_min_corr(). */
 double reserve_reward( Index t , double p , double H ) const;

 /// the on-to-on transition (4) into instant t, reserves included
 /** Computes in \p out the function \f$ G \f$ of (4) on
  * \f$ [ \mbox{lo} , \mbox{hi} ] \f$, where the ramp terms of the band
  * (3) are those of \p ramp_up and \p ramp_down and the window of the move
  * is \f$ [ -\omega^- , \omega^+ ] \f$, \f$ \omega^+ \f$ being \p win_up and
  * \f$ \omega^- \f$ being \p win_down. When no reserve cost is negative at
  * \p t this is the plain sliding minimum [see sliding_min()].
  *
  * Otherwise \f$ G \f$ is built exactly, piece by piece, on the following
  * grounds. Let \f$ \bar{\rho}_1 \f$ be the fraction \f$ \rho \f$ of the
  * reserve with the most negative cost and, if both costs are negative,
  * \f$ \bar{\rho}_2 = \rho^{pr}_t + \rho^{sc}_t \f$, so that
  * \f$ H \mapsto g_t( p , H ) \f$ in (2) is piecewise linear, with the
  * most negative cost as its slope on \f$ [ 0 , \bar{\rho}_1 p ] \f$, the
  * other cost on \f$ [ \bar{\rho}_1 p , \bar{\rho}_2 p ] \f$ (if it is
  * negative) and 0 beyond. For a
  * fixed \f$ p \f$ the function
  * \f$ \phi_p( q ) = F( q ) + g_t( p , \mathcal{H}_t( p , q ) ) \f$ is convex
  * on the window (see the file comment), and on every set of \f$ p \f$ where
  * the configuration of a minimizer \f$ q^*( p ) \f$ is fixed (the piece of
  * \f$ F \f$ holding it, the term of (3) that gives \f$ H \f$, the segment
  * of \f$ g_t \f$ and the active end of the window, if any) the optimality
  * condition \f$ 0 \in \partial \phi_p( q ) \f$ is affine in
  * \f$ ( p , q ) \f$. Hence \f$ q^* \f$ is either constant there (a kink of
  * \f$ F \f$, the stationary point of a quadratic piece of \f$ F \f$ plus a
  * slope of \f$ g_t \f$, an end of the domain of \f$ F \f$, or one of
  * \f$ \Delta^-_t + P^{mn}_t \f$ and \f$ K_t - \Delta^+_t \f$, where a ramp
  * term of (3) equals the band) or it runs on one of the lines
  * \f[
  *   q = p - \omega^+ \, , \;\; q = p + \omega^- \, , \;\;
  *   q = p - \tfrac{1}{2} ( \Delta^+_t - \Delta^-_t ) \, , \;\;
  *   q = 2 p - \Delta^+_t - P^{mn}_t \, , \;\;
  *   q = 2 p + \Delta^-_t - K_t \, , \;\;
  *   q = ( 1 \pm \bar{\rho}_m ) p \mp \Delta^\pm_t \tag{5}
  * \f]
  * (an end of the window, the peak of the two ramp terms of (3), a ramp
  * term equal to the band, a ramp term equal to \f$ \bar{\rho}_m p \f$,
  * \f$ m = 1 , 2 \f$, with \f$ \Delta^+_t \f$ for the upper sign and
  * \f$ \Delta^-_t \f$ for the lower one). On such a set
  * \f$ q^*( p ) = s p + c \f$, \f$ F \f$ is one quadratic and \f$ g_t \f$ is
  * affine in \f$ ( p , H ) \f$ with \f$ H \f$ affine in \f$ ( p , q ) \f$,
  * so that \f$ G( p ) = \phi_p( s p + c ) \f$ is a quadratic in \f$ p \f$
  * with curvature \f$ s^2 \f$ times that of the piece of \f$ F \f$; and
  * \f$ G \f$ is convex, being the minimum over \f$ q \f$ of a function
  * jointly convex in \f$ ( p , q ) \f$ on the convex set
  * \f$ p - \omega^+ \leq q \leq p + \omega^- \f$. The construction sweeps
  * \f$ p \f$ upwards: it computes \f$ q^* \f$ exactly at two points just right
  * of the current \f$ p \f$ and takes as the current set the line of (5)
  * through both points, or the constant \f$ q^* \f$ if the measured slope is
  * below \f$ 1 / 2 \f$ (a measured locus that matches no line is kept as it
  * is). Then it emits the quadratic piece up to the next point where the
  * configuration can change, and accepts it only if the optimality
  * conditions of \f$ q^* \f$ hold at both ends of it; since these
  * conditions are affine in \f$ p \f$ under a fixed configuration, they
  * then hold in between. A slope
  * below \f$ 1 / 2 \f$ is that of a constant only when every
  * \f$ \bar{\rho}_m \f$ is at most \f$ 1 / 2 \f$: the line
  * \f$ q = ( 1 - \bar{\rho}_m ) p + \Delta^-_t \f$ of a larger
  * \f$ \bar{\rho}_m \in ( 1 / 2 , 1 ) \f$ (e.g., with both reserves
  * rewarded and the ramp-down term of (3) binding) has a slope in
  * \f$ ( 0 , 1 / 2 ) \f$, and it is taken whenever it fits the two points
  * within the relative tolerance \f$ 10^{-6} \f$ and better than the
  * constant does. The lines of a \f$ \bar{\rho}_m \geq 1 \f$, whose slope
  * is not positive, are never tried, and need not be: if \f$ P^{mn}_t > 0
  * \f$ such a locus is never active, since \f$ \mathcal{H}_t \leq p -
  * P^{mn}_t < \bar{\rho}_m p \f$, i.e., the fraction never binds before the
  * band, and if \f$ P^{mn}_t = 0 \f$ and \f$ \bar{\rho}_m = 1 \f$ the line
  * \f$ q = \Delta^-_t \f$ is the constant \f$ \Delta^-_t + P^{mn}_t \f$
  * above, and the other one the line of the band term. Wherever no set is
  * identified, or the budget of steps of the sweep is spent, the rest of
  * the domain is covered by linear pieces that interpolate exact values of
  * \f$ G \f$, refined until the value at the midpoint is within
  * \f$ 10^{-7} \max\{ 1 , | G | \} \f$ of the secant (\f$ G \f$ at the
  * midpoint of the part covered) or the piece is shorter than
  * \f$ 10^{-5} \f$, within a budget of \f$ 2 \cdot 10^5 \f$ bisections that
  * the cases met so far never reach, so that the domain of \f$ G \f$, i.e.,
  * the set of reachable powers at \p t, is not cut short, and \f$ G \f$ is
  * approximated there from above (a secant of a convex function).
  *
  * @param F the value function at \f$ t - 1 \f$, convex
  * @param ramp_up \f$ \Delta^+_t \f$ in the ramp terms of (3)
  * @param ramp_down \f$ \Delta^-_t \f$ in the ramp terms of (3)
  * @param lo the lower end of the domain of \f$ G \f$
  * @param hi the upper end of the domain of \f$ G \f$
  * @param t the instant of the landing power
  * @param out the function \f$ G \f$ (cleared first)
  * @param acap the cap \f$ K_t \f$ of (3): \f$ P^{mx}_t \f$ if negative
  *        (an instant inside an on-run), the smaller of \f$ P^{mx}_t \f$
  *        and \f$ P^{sd}_{t+1} \f$ when the run closes at \p t
  * @param win_up \f$ \omega^+ \f$, \p ramp_up if NaN (the window of a thermal
  *        unit)
  * @param win_down \f$ \omega^- \f$, \p ramp_down if NaN; a derived unit whose
  *        scheduled move is narrower than the ramp passes a smaller window,
  *        the ramp terms of (3) staying those of the ramp */
 void sliding_min_corr( const PQFun & F ,
                        double ramp_up , double ramp_down ,
                        double lo , double hi , Index t , PQFun & out ,
                        double acap = -1.0 ,
                        double win_up = std::nan( "" ) ,
                        double win_down = std::nan( "" ) );

 /// the minimizer \f$ q \f$ of (4) at a fixed power \p p at t
 /** Returns a \f$ q \f$ that attains the minimum in (4) at \p p, with the
  * same arguments as sliding_min_corr(), i.e., the power at \f$ t - 1 \f$
  * of the schedule the value \f$ G( p ) \f$ priced; the backward pass of
  * the solvers uses it to recover the schedule. */
 double reserve_corr_argmin( const PQFun & F , double ramp_up ,
                             double ramp_down , Index t , double p ,
                             double acap = -1.0 ,
                             double win_up = std::nan( "" ) ,
                             double win_down = std::nan( "" ) ) const;

/*--------------------------------------------------------------------------*/
/*---------------- PIECEWISE-QUADRATIC FUNCTION HELPERS --------------------*/
/*--------------------------------------------------------------------------*/

 /// evaluate a quadratic piece at p
 static double eval_piece( const PieceQuad & pc , double p ) {
  return( pc.alfa * p * p + pc.beta * p + pc.gamma );
 }

 /// evaluate a PQFun F at p; TUEDPINF outside the domain of F
 static double eval( const PQFun & F , double p );

 /// minimizer of a single quadratic piece on its interval
 static double argmin_piece( const PieceQuad & pc );

 /// minimum of a PQFun over [ lo , hi ] intersected with its domain
 /** Returns the pair (minimum, minimizer), and ( TUEDPINF , 0 ) if the
  * intersection is empty. */
 static std::pair< double , double > min_over(
  const PQFun & F , double lo , double hi );

 /// pointwise addition of the constant \p c to F
 static void shift_by( PQFun & F , double c );

 /// pointwise addition of \f$ \alpha p^2 + \beta p + \gamma \f$ to F
 static void add_quadratic( PQFun & F ,
                            double alfa , double beta , double gamma );

 /// restrict F to the domain [ lo , hi ] (in place)
 static void clamp_domain( PQFun & F , double lo , double hi );

 /// add a piecewise linear PQFun G to F on the domain of F, in place
 /** Pointwise sum \f$ F( p ) \leftarrow F( p ) + G( p ) \f$ on the domain
  * of F, whose pieces are split at the breakpoints of G; outside the domain
  * of G, G counts as 0. */
 static void add_pwq( PQFun & F , const PQFun & G );

 /// true if F1 is pointwise not smaller than F2, up to \p eps, on the whole
 /// domain of F1; false if the domain of F1 is not within that of F2
 static bool is_dominated_by( const PQFun & F1 , const PQFun & F2 ,
                              double eps = 1e-9 );

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

 /// the sliding minimum of a convex PQFun F over the window of a move
 /** Computes in \p out (cleared first) the function
  * \f[
  *   G( p ) = \min \{ \, F( q ) \, : \, p - \omega^+ \leq q \leq p + \omega^-
  *   \, \}
  * \f]
  * on \f$ [ \mbox{lo} , \mbox{hi} ] \f$, \f$ \omega^+ \f$ being \p ramp_up
  * and \f$ \omega^- \f$ being \p ramp_down, which is (4) with no reserve. With
  * \f$ p^* \f$ a minimizer of F, \f$ G \f$ is F shifted rightwards by
  * \f$ \omega^+ \f$ on the right of \f$ p^* \f$, leftwards by \f$ \omega^- \f$
  * on its left, and constant at \f$ F( p^* ) \f$ in between; a window that
  * does not contain 0 (a move that must go upwards, or downwards) is
  * reduced to one that does by a translation, and the window of a single
  * point is the translation of F. It uses a scratch buffer of the object,
  * hence it is not static. */
 void sliding_min( const PQFun & F ,
                   double ramp_up , double ramp_down ,
                   double lo , double hi , PQFun & out );

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
 std::vector< double > shutdown_costs;   ///< empty if shutting down is free
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 // whether the Block actually defines ramp limits (as opposed to defaulting
 // them to the largest of InitialPower and MaxPower, see
 // load_common_parameters()); the shut-down at 0 of a unit on before
 // the horizon is allowed iff InitialPower <= ShutDownLimit[ 0 ], whatever
 // the ramps, in both solvers as in every formulation
 bool has_ramp_up{ false };
 bool has_ramp_down{ false };
 std::vector< double > min_power;
 std::vector< double > max_power;
 std::vector< double > bound_on;
 std::vector< double > bound_down;

 std::vector< double > quad_term;
 std::vector< double > linear_term;
 std::vector< double > const_term;

 // -- reactive power --------------------------------------------------- //
 std::vector< double > reactive_linear_term;
 std::vector< double > reactive_min;
 std::vector< double > reactive_max;
 std::vector< double > reactive_min_on;
 std::vector< double > reactive_max_on;

 // -- spinning reserve ------------------------------------------------- //
 // participation factors (caps in pr<=rho_p*p, sr<=rho_s*p) and objective
 // cost coefficients on the reserve variables; each empty if the
 // corresponding reserve is absent. The cost may be a Lagrangian price
 // (possibly negative, i.e. a reward), kept separate from the rho cap.
 std::vector< double > primary_rho;
 std::vector< double > secondary_rho;
 std::vector< double > primary_reserve_cost;
 std::vector< double > secondary_reserve_cost;

 // -- design (investment) --------------------------------------------- //
 bool   has_design{ false };  ///< true iff the unit has an investment cost
 double design_cost{ 0 };     ///< design variable objective coefficient
 bool   design_on{ false };   ///< the design decision computed by run_DP()

 // -- fixed Variable handling ----------------------------------------- //
 // the commitment (and design) Variable can be fixed, which run_DP() honors
 // by killing the ON states of the instants fixed OFF and the OFF states of
 // the instants fixed ON (see load_fixings()); this may make the problem
 // infeasible.
 std::vector< Index > nxt_off;  ///< first instant >= t fixed OFF (T if none)
 std::vector< Index > nxt_on;   ///< first instant >= t fixed ON (T if none)
 bool f_has_fixings{ false };   ///< true iff some commitment is fixed
 bool fixed_to_max{ false };    ///< "FixToMaximum" of the unit is positive
 bool f_must_build{ false };    ///< commitment fixed ON, or design fixed to 1
 bool f_no_build{ false };      ///< the design variable is fixed to 0

 double eps{ 1e-10 };         ///< numerical tolerance

 // -- reusable scratch for the PQ machinery --------------------------- //
 // static thread_local so the per-source sweeps of the parallel base
 // solver (parDP) can call sliding_min()/sliding_min_corr() concurrently
 // without racing on shared scratch; for the single-threaded run-length
 // solver the semantics are identical (one thread, storage merely recycled)
 static thread_local std::vector< PQFun > m_pqpool;  ///< spare PQFun storage
 static thread_local PQFun                m_raw;      ///< sliding_min scratch

/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitDPSolverBase ) )

/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*----------------- diagnostic profiling counters --------------------------*/
/*--------------------------------------------------------------------------*/
// Only compiled when TUEDPS_PROFILE > 0 (dev instrumentation, 0 in
// production). These counters are shared by the reserve machinery /
// sliding_min_corr() (ThermalUnitDPSolverBase.cpp) and compute()'s
// cumulative PROF print (ThermalUnitExtDPSolver.cpp), so they need a single
// storage across the two translation units; inline variables give exactly
// one definition. Each of those two .cpp #defines TUEDPS_PROFILE *before*
// including this header, so this block is compiled only there, and only
// when profiling is on.

#if defined( TUEDPS_PROFILE ) && TUEDPS_PROFILE
namespace SMSpp_di_unipi_it
{
 inline unsigned long g_smc = 0 , g_gmin = 0 , g_node = 0 ,
                      g_Fsum = 0 , g_Fmax = 0 , g_pieces = 0 ,
                      g_gmin_bite = 0 , g_smc_bite = 0 , g_fast = 0 ,
                      g_param_calls = 0 , g_param_bad = 0 , g_param_pcs = 0;
 inline double g_param_maxdev = 0 , g_param_ormaxdev = 0 , g_gmin_maxdev = 0;
}
#endif

#endif  /* ThermalUnitDPSolverBase.h included */

/*--------------------------------------------------------------------------*/
/*--------------- End File ThermalUnitDPSolverBase.h -----------------------*/
/*--------------------------------------------------------------------------*/
