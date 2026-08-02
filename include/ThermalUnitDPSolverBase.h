/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitDPSolverBase.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolverBase class, the shared base of the
 * two dynamic-programming Solvers for the single-Unit Commitment (1UC)
 * problem on a ThermalUnitBlock: ThermalUnitDPSolver (the on/off-graph +
 * per-source economic-dispatch sweep) and ThermalUnitExtDPSolver (the
 * run-length "multi-layer ON / single-layer OFF" scheme).
 *
 * The two solvers differ only in how they structure the dynamic program;
 * they share, and this base class owns:
 *
 *  - all the data loaded from the ThermalUnitBlock (power/ramp/cost bounds,
 *    reactive box, spinning-reserve participation factors and prices, design
 *    and commitment-fixing data), plus the common loader
 *    load_common_parameters();
 *
 *  - the convex piecewise-quadratic value-function type (PieceQuad / PQFun)
 *    and every "shape" operation on it (eval, argmin, sliding minimum,
 *    pointwise add, clamp, dominance test, storage pooling);
 *
 *  - the whole spinning-reserve model: the per-period greedy reserve LP
 *    (reserve_alloc / reserve_alloc_band / reserve_reward), the convex
 *    piecewise-linear reserve discount \f$ g_t \f$
 *    (build_reserve_discount), the residual-ramp on->on transition
 *    sliding_min_corr() (which folds the deliverability-correct reserve
 *    reward into the ramp-window minimisation), and the corr-aware
 *    predecessor recovery reserve_corr_argmin().
 *
 * None of this machinery touches either solver's DP tables: the value
 * function is always passed in as a PQFun and every quantity read is one of
 * the loaded vectors below, so the two solvers can share a single
 * implementation.
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

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ThermalUnitDPSolverBase -----------------------*/
/*--------------------------------------------------------------------------*/
/// shared base of the two 1UC dynamic-programming Solvers
/** Holds the data loaded from the ThermalUnitBlock, the convex
 * piecewise-quadratic value-function machinery and the whole spinning-reserve
 * model (including the residual-ramp on->on transition). It derives from
 * Solver but is abstract: the concrete DP structure, and the Solver interface
 * (compute(), get_var_solution(), ...), are provided by the two subclasses
 * ThermalUnitDPSolver and ThermalUnitExtDPSolver. */

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
/*------------------------- PROTECTED TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// one quadratic piece \f$ \alpha p^2 + \beta p + \gamma \f$ on the
 /// interval [left, right] of the power variable \f$ p \f$
 /** The atomic building block of the piecewise-quadratic value functions
  * used by the DP. The three coefficients encode the full univariate
  * quadratic \f$ \alpha p^2 + \beta p + \gamma \f$; the two endpoints
  * [left, right] delimit the portion of the power axis on which this
  * expression is meaningful. Neighbouring pieces in a PQFun share their
  * adjacent endpoints (right of piece i == left of piece i+1). */

 struct PieceQuad {
  double alfa;   ///< coefficient of \f$ p^2 \f$ (>= 0 for a convex piece)
  double beta;   ///< coefficient of \f$ p \f$
  double gamma;  ///< additive constant (absorbs accumulated path cost)
  double left;   ///< left endpoint of the piece (inclusive)
  double right;  ///< right endpoint of the piece (inclusive)
 };

 /// a convex piecewise quadratic function on the power axis
 /** Stored as a vector of PieceQuad sorted by left endpoint and covering a
  * contiguous sub-interval of the power axis with no gaps and no overlaps
  * (except at shared endpoints). The invariants are:
  *
  *  1. pieces[ i ].right == pieces[ i+1 ].left for all i (continuous
  *     support);
  *  2. alfa >= 0 on every piece (each piece is itself convex);
  *  3. the values at shared endpoints agree (function is continuous);
  *  4. the subgradient is monotonically non-decreasing across endpoints,
  *     so the *overall* piecewise function is convex.
  *
  * The empty vector represents the +INFinity function, used to signal an
  * infeasible or not-yet-reached DP state. All the "shape" operations
  * preserve invariants (1)-(4) by construction. */

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

 /// expand a per-instant vector: a single value is broadcast over the horizon
 void retrieve_term( std::vector< double > & out ,
                     const std::vector< double > & in ) const;

/*--------------------------------------------------------------------------*/
/*---------------------------- RESERVE MODEL -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// optimal spinning-reserve provision at instant t given active power p
 /** Solves, for the unit on at \p t with power \f$ p \f$, the per-period
  * reserve LP
  * \f[
  *  \min \{ c_{pr} \, pr + c_{sr} \, sr \, : \,
  *          pr \leq \rho_p \, p \, , \, sr \leq \rho_s \, p \, , \,
  *          pr + sr \leq cap - p \, , \, pr \, , \, sr \geq 0 \}
  * \f]
  * (reserve only helps when its cost coefficient is negative, i.e. a
  * Lagrangian reward). Writes the optimal \f$ pr \f$, \f$ sr \f$ and
  * returns the optimal value \f$ g_t( p ) = c_{pr} \, pr + c_{sr} \, sr \f$
  * (0 in the standalone cost case). */
 double reserve_alloc( Index t , double p , double & pr , double & sr ,
                       double cap ) const;

 /// like reserve_alloc() but with the shared reserve band H passed directly
 /** Greedily fills the reserves against the band \p H (rather than deriving
  * it from a power cap), writing the optimal pr, sr and returning the
  * reward. */
 double reserve_alloc_band( Index t , double p , double H ,
                            double & pr , double & sr ) const;

 /// reserve "discount" \f$ g_t( p ) \f$ as a convex piecewise-linear PQFun
 /** Builds \f$ g_t( p ) \f$ = min over \f$ ( pr , sr ) \f$ in the reserve
  * polytope of \f$ c_{pr} \, pr + c_{sr} \, sr \f$, as a function of the
  * production \f$ p \f$, on [min_power[t], cap]. Non-positive (a reward)
  * when some reserve price is negative, and identically zero (empty PQFun)
  * when no price is negative.
  * \p cap is the upper power cap of the reserve band: max_power[t] at an
  * interior period, the tighter start-up cap bound_on[t] or shut-down cap
  * bound_down[t] at a boundary period. */
 PQFun build_reserve_discount( Index t , double cap ) const;

 /// greedy spinning-reserve reward at power \p p under a shared band \p H
 /** The optimal (most negative) value of the per-period reserve LP when the
  * reserve room (the shared band) is \p H. Like reserve_alloc() but with the
  * band passed directly rather than derived from a power cap; used by
  * sliding_min_corr(). Returns 0 when \p H <= 0 or no reserve is priced. */
 double reserve_reward( Index t , double p , double H ) const;

 /// on->on transition with the residual-ramp reserve penalty folded in
 /** Computes
  * \f[
  *  out( p ) = \min \{ F( q ) + corr_t( q , p ) \, : \,
  *                     q \in [ p - ramp\_up , p + ramp\_down ] \}
  * \f]
  * on \f$ p \in [ lo , hi ] \f$, where \f$ corr_t( q , p ) \f$ is the
  * reserve reward under the residual-ramp band \f$ \min( A_t( p ) ,
  * B( p - q ) ) \f$, with \f$ A_t( p ) = \min( p - min\_power , cap - p )
  * \f$ the capacity band and \f$ B( d ) = \min( ramp\_up - d ,
  * ramp\_down + d ) \f$ the ramp tent. When no reserve is rewarded
  * corr == 0 and this falls back to the exact sliding_min().
  * @param acap upper cap of the capacity band \f$ A_t \f$: defaults
  * (acap < 0) to max_power[t] at an interior step, set to the shut-down cap
  * bound_down[t+1] when the transition closes a run. */
 void sliding_min_corr( const PQFun & F ,
                        double ramp_up , double ramp_down ,
                        double lo , double hi , Index t , PQFun & out ,
                        double acap = -1.0 );

 /// argmin over the ramp window of \f$ F( q ) + corr_t( q , p ) \f$ at a
 /// fixed landing \f$ p \f$
 /** The corr-aware predecessor power for the on->on transition into \p t:
  * the \f$ q \f$ that minimises \f$ F( q ) + corr_t( q , p ) \f$ over
  * \f$ [ p - ramp\_up , p + ramp\_down ] \cap dom( F ) \f$. Used by the
  * backward pass so the recovered power profile is the one the DP value
  * priced. Returns the energy argmin when no reserve is rewarded at \p t. */
 double reserve_corr_argmin( const PQFun & F , double ramp_up ,
                             double ramp_down , Index t , double p ,
                             double acap = -1.0 ) const;

/*--------------------------------------------------------------------------*/
/*---------------- PIECEWISE-QUADRATIC FUNCTION HELPERS --------------------*/
/*--------------------------------------------------------------------------*/

 /// evaluate a quadratic piece (alfa, beta, gamma) at p
 static double eval_piece( const PieceQuad & pc , double p ) {
  return( pc.alfa * p * p + pc.beta * p + pc.gamma );
 }

 /// evaluate a piecewise quadratic function F at p; TUEDPINF outside dom(F)
 static double eval( const PQFun & F , double p );

 /// argmin of a single quadratic piece on [left, right]
 static double argmin_piece( const PieceQuad & pc );

 /// minimum of a PQFun over [lo, hi] intersected with its domain
 /** Returns { min_value, argmin_p }; { TUEDPINF, 0 } if the intersection is
  * empty. */
 static std::pair< double , double > min_over(
  const PQFun & F , double lo , double hi );

 /// pointwise addition of the constant 'c' to F
 static void shift_by( PQFun & F , double c );

 /// pointwise addition of \f$ \alpha p^2 + \beta p + \gamma \f$ to F
 static void add_quadratic( PQFun & F ,
                            double alfa , double beta , double gamma );

 /// restrict F to the domain [lo, hi] (in place)
 static void clamp_domain( PQFun & F , double lo , double hi );

 /// add a (piecewise-linear) PQFun G to F on F's domain, in place
 /** Pointwise sum \f$ F( p ) \leftarrow F( p ) + G( p ) \f$ for
  * \f$ p \in dom( F ) \f$; pieces of F are split where a breakpoint of G
  * falls. Where G is undefined it contributes 0. */
 static void add_pwq( PQFun & F , const PQFun & G );

 /// check whether F1 is pointwise >= F2 (up to tolerance eps) on all of
 /// dom( F1 ); false if dom( F1 ) extends beyond dom( F2 ).
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

 /// ramp-constrained sliding minimum of a convex piecewise quadratic F
 /** Computes
  * \f[
  *  G( p_t ) = \min \{ F( q ) \, : \,
  *                     q \in [ p_t - ramp\_up , p_t + ramp\_down ] \}
  * \f]
  * on the domain [lo, hi]. Assumes F convex. The result is written into
  * \p out (cleared first); an internal scratch buffer (m_raw) is reused
  * across calls, so the hot path performs no per-call allocation.
  * Non-static for that reason. */
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
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 // whether the Block actually defines ramp limits (as opposed to defaulting
 // them to max_power); the start-up / shut-down trajectory of the initial
 // state is only enforced when the corresponding ramp limit is present
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
