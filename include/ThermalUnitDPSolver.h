/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolver class, that solves the
 * ThermalUnitBlock, comprised its spinning-reserve variables, using a
 * Dynamic Programming algorithm.
 *
 * \author Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Claudio Gentile, Antonio Frangioni, Niccolo' Iardella,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ThermalUnitDPSolver
 #define __ThermalUnitDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <memory>

#include "Solver.h"

#include "ThermalUnitBlock.h"

#include "ThermalUnitDPSolverBase.h"

#define TUDPS_PARALLEL 1
/* If TUDPS_PARALLEL > 0, the (independent) per-ON-node Economic Dispatch
 * solves in compute_EDPs() can be run in parallel with FastFlow; the actual
 * number of workers is controlled at run time by the intMaxThread parameter
 * (0 = use all available cores, 1 = serial, k = k workers). Set to 0 in build
 * setups where FastFlow is not available (e.g. the plain makefiles). */

#ifndef TUDPS_PAR_MIN_N
 #define TUDPS_PAR_MIN_N 768
#endif
/* Default time horizon below which compute_EDPs() stays serial even with
 * intMaxThread != 1: the thread-dispatch overhead is not amortised on short
 * horizons (benchmarks on 8 cores put the serial/parallel break-even around
 * 600 time steps, so this is set conservatively above it). It is only the
 * default of the run-time intParMinN parameter (see below); set intParMinN in
 * the ComputeConfig to override it per solver instance (0 forces the parallel
 * path on every horizon), or override this compile-time default with
 * -DTUDPS_PAR_MIN_N=<n>. */

#if TUDPS_PARALLEL
namespace ff { class ParallelFor; }  // forward declaration (FastFlow)
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitDPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving a Single Unit Commitment problem with a DP approach
/** The ThermalUnitDPSolver is a Solver for tackling the single-Unit
 * Commitment problem with ramp-up and ramp-down (as well as minimum up-
 * and down-time) constraints and convex quadratic separable objective.
 * The solver uses a Dynamic Programming approach, recasting the problem as
 * a shortest path one on a graph with the following structure
 *
 *       (0,1)  (1,1)  (2,1)  (3,1)  (4,1)  (5,1)  ...  (n-1,1)
 *  (s)                                                           (d)
 *       (0,0)  (1,0)  (2,0)  (3,0)  (4,0)  (5,0)  ...  (n-1,0)
 *
 * n being the length of the time horizon.
 *
 * The fundamental tool for solving this problem is the ability of efficiently
 * solving Economic Dispatch problems that find the min-cost energy production
 * of the unit if it is on for a continuous time interval. In particular we
 * denote by ED( h , k ) the total power cost (but not the fixed and start-up
 * ones, that are computed separately) of the unit if started up exactly at
 * the beginning of time h >= 0 and shut down exactly at the end of time
 * h <= k <= n - 1, i.e., being online for all the time instants h, h + 1,
 * ..., k (note that h = k is possible). Similarly, we denote by SUC( h , k )
 * the cost of having the unit off from the beginning of h to the end of k
 * and then starting up at k + 1: this is typically easy to compute.
 *
 * The problem is therefore reduced to a Shortest Path between (s) and (d)
 * on the acyclic graph constructed as follows:
 *
 * - Each arc from an ON node ( i , 1 ) to an OFF node ( j , 0 ), for
 *   0 <= i < j <= n - 1, means that the unit is started at the beginning of
 *   time i and shut down at the end of time j - 1, so that it is off at
 *   time j. The cost of the arc is therefore ED( i , j - 1 ). Note that
 *   this problem concerns the power variables p[ i ], p[ i + 1 ], ...
 *   p[ j - 1 ], but also implicitly p[ j ] that will be necessarily
 *   fixed to 0 (as the unit is down at j). However, such an arc exists only
 *   if j - i = number of consecutive periods the unit remains on is
 *   >= min up-time. In particular, j == i + 1, i.e., the unit is started
 *   up and shut down at the end of the same period, is only possible if
 *   min up-time <= 1, i.e., there is no min up-time requirement. Note that
 *   min up-time need necessarily be >= 1 as the unit cannot remain on for
 *   less than one time instant, so min up-time == 0 hardly makes sense.
 *
 * - Each arc from an OFF node ( i , 0 ) to an ON node ( j , 1 ), for
 *   0 <= i < j <= n - 1, means that the unit is shut down at the beginning
 *   of time i and remains down up until the end of time j - 1, then it is
 *   started up at j. Hence, the cost of the arc is the (possibly,
 *   time-variable) start-up costs SUC( i , j - 1 ). This basically fixes
 *   p[ i ] = p[ i + 1 ] = ... = p[ j - 1 ] = 0, but leaves p[ j ] free to be
 *   anything (it will be decided by the outgoing arcs of ( j , 1 )).
 *   Such an arc exists only if j - i = number of consecutive periods the
 *   unit remains off is >= min down-time. In particular, in this case it
 *   would even be possible i == j, i.e., the unit was shut down at the
 *   end of period i - 1 and it immediately re-started at the beginning
 *   of period i, if min down-time == 0, i.e., there is no min down-time
 *   requirement. Note that, unlike for min up-time, in this case the
 *   value 0 in principle makes sense and it is different from the value
 *   1. However, we avoid "vertical" arcs ( i , 0 ) --> ( i , 1 ) since
 *   a solution where the unit is shut down and immediately started up is
 *   never economical w.r.t. one where the unit is never shut down in the
 *   first place, since shutting down entails a "shutdown trajectory" that
 *   brings the unit to the right stopping power that further constrains
 *   the unit (but this is not prohibited if the unit remains on, which
 *   means that remaining on is always at least as cheap).
 *
 * - From each ON node ( i , 1 ) there always is one arc to the destination
 *   d, meaning that the unit remains on in all the time instants between i
 *   and n - 1, and it is *not* shut down at the end of the period. The
 *   cost of this arc is the optimal cost of a "special" ED( i , n - 1 ),
 *   deciding on all variables p[ i ], p[ i + 1 ], ..., p[ n - 1 ] and
 *   *not* (implicitly) fixing p[ n - 1 ] = 0 as ED( i , n - 2 ),
 *   corresponding to the arc ( i , 1 ) --> ( n - 1 , 0 ) does. The reason
 *   why ED( i , n - 1 ) is "special" is that, due to the ramp-down
 *   constraints, if the unit has to be down at time k, then it must enter
 *   in a "shutdown trajectory" in the previous time instants, so that the
 *   final power p[ k - 1 ] is the right one to stop. This constrains the ED,
 *   resulting in a higher cost. This means that forcing the shut down at the
 *   end of n - 1 is never economical: it is in principle better to allow the
 *   unit do what it wants (which may comprise autonomously entering in a
 *   shutdown trajectory if this is the optimal thing to do, as this is not
 *   prohibited). This is why the constraints of ED( h , n - 1 ) do not
 *   include the one forcing the power of the unit at the last time instant to
 *   be the shutdown one, unlike for all the other ED( h , k ).
 *
 * - From each OFF node ( i , 0 ) there always is one arc to the destination
 *   d, meaning that the unit remains off in all the time instants between i
 *   and n - 1. Ordinarily this would imply that the unit is started up right
 *   at the beginning of the next horizon of operations, but this is not of
 *   our concern for the current problem. This means that all these arcs have
 *   *zero cost*, as any startup cost will be accounted for in the next
 *   horizon of operations, if any.
 *
 * - From the node (s) there are arcs going to either ( i , 1 ) or ( i, 0 )
 *   nodes depending on the value of init_up_down_time (the amount of time
 *   the unit has been on or off prior to the initial time instant 0), the
 *   min_up_time, min_down_time, initial_power and delta_ramp_down values,
 *   as applicable. The rules are the following:
 *
 *   = If init_up_down_time > 0, then the unit has been on for
 *     init_up_down_time periods before the initial time instant 0 and it
 *     is at power initial_power at the beginning of time instant 0. Hence,
 *     by the ramp-down constraints, there is a minimum number of time
 *     instants, t_ramp_min, that are necessary to bring the unit to the
 *     power level required to stop (t_ramp_min could be 0 if initial_power
 *     happens to be exactly the right power level). Then, there will be
 *     arcs between s and all nodes ( i , 0 ) for "sufficiently large"
 *     i >= min_node, where:
 *
 *     * if init_up_down_time >= min_up_time, i.e., the unit is already
 *       on since long enough to satisfy the min up-time constraint, then
 *       min_node = t_ramp_min;
 *
 *     * if init_up_down_time < min_up_time, i.e., the unit has to remain
 *       on anyway for at least other ( min_up_time - init_up_down_time )
 *       instants, then
 *       min_node = max( t_ramp_min , min_up_time - init_up_down_time );
 *
 *     The cost of each arc ( s , i ) will be ED( 0 , i - 1 ), where these
 *     ED are also "special" in the sense that, for the sake of ramp-up and
 *     ramp-down constraints, the initial_power value is used as reference.
 *     Note however the very special case where init_up_down_time >=
 *     min_up_time and t_ramp_min == 0, i.e., min_node == 0. This
 *     corresponds to the case where the unit is "on but on the brink of
 *     shutting down" at the beginning of the time horizon. This means that
 *     there will be *both* an arc s --> ( 0 , 1 ) saying "the unit is on
 *     and will remain on for a while", *and* an arc s --> ( 0 , 0 ) saying
 *     "the unit is on but I'll shut it down immediately and will remain
 *     off for a while".
 *
 *     Note that "sufficiently large" i includes i == n, i.e., the arc
 *     s --> d corresponding to "the unit was on at the beginning and
 *     remains on for the whole period".
 *
 *   = If init_up_down_time <= 0, then the unit has been off for
 *     init_up_down_time periods before the initial time instant 0; note
 *     that the value 0 is included, meaning "the unit has just been
 *     shut down when the time begins". Then, there will be arcs between
 *     s and all nodes ( i , 1 ) for "sufficiently large"
 *     i >= max( min down-time + init_up_down_time , 0 ). That is, we
 *     wait for the remaining min down-time - ( - init_up_down_time )
 *     time periods (if any) required by the min down-time constraint
 *     before allowing the unit to be started up again. These arcs will
 *     have cost SUC( 0 , i ) since the unit will be started up at i.
 *     Note that if min down-time == 0, i.e., there is no min
 *     down-time requirement, this means that i == 0 is always possible,
 *     i.e., the unit is started up immediately at the beginning. This
 *     potentially yields the "double strange" case where
 *     init_up_down_time == 0, i.e., the unit had just been shut down
 *     and it is immediately restarted. While this is not economical, we
 *     cannot (and have no reason to) avoid it, as in the case of
 *     OFF -> ON arcs, because the decision to shut down the unit right
 *     at the end of the previous interval (encoded by init_up_down_time
 *     == 0) is not in our hands as it has been taken before our time has
 *     come. Yet, this is neither impossible nor logically contradictory, so
 *     there is no problem (and min down-time == 0 is unlikely anyway).
 *
 *     Note that "sufficiently large" i includes i == n, i.e., the arc
 *     s --> d corresponding to "the unit was off at the beginning and
 *     remains off for the whole period". Like all arcs OFF -> d, this
 *     does not really imply that the unit will necessarily be restarted
 *     immediately after, and anyway this is outside the boundaries of
 *     the current problem; hence, this arc has *zero cost*.
 *
 * ThermalUnitDPSolver first builds the graph, then uses one EDSolver for
 * each ON node (comprised s if the unit is on at the beginning, and
 * therefore it is equivalent to an ON node) to solve EDs to compute the arc
 * costs, then uses an (acyclic) min-path algorithm to solve the problem. */

class ThermalUnitDPSolver : public ThermalUnitDPSolverBase
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 static constexpr auto TUDPINF = Inf< double >();  ///< the INF value

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 * @{ */

 ThermalUnitDPSolver( void ) {};

 ~ThermalUnitDPSolver() override;  // defined in the .cpp (pimpl'd FastFlow)

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public methods derived from base classes
 * @{ */

 /// sets the Block that the Solver has to solve
 void set_Block( Block * block ) override;

 /// solves the constructed problem
 int compute( bool changedvars = true ) override;

 /// tells whether a solution is available
 bool has_var_solution( void ) override { return( f_end.pred ); }

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc ) override;

/*--------------------------------------------------------------------------*/
 /// returns the schedule the DP has found as a ThermalUnitBlockSolution
 /** Returns the schedule the dynamic programming has found as a
  * ThermalUnitBlockSolution [see ThermalUnitBlock.h], filled straight out of
  * the data structures of the Solver rather than by writing it in the
  * Variable of the ThermalUnitBlock and having it read back from there: no
  * abstract representation is therefore required to exist, and the
  * ThermalUnitBlock is not written into at all, hence it is not lock()-ed
  * and any number of Solver attached to it can produce their own Solution at
  * the same time.
  *
  * The active power, the commitment and the dimensioning variable are always
  * saved; the spinning reserves and the reactive power only if the unit has
  * them, which is the same rule get_var_solution() follows by only writing
  * the Variable that exist. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// the Solution is filled from the data of the DP, not from the Variable

 [[nodiscard]] bool is_get_Solution_physical( void ) const override {
  return( true );
  }

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override { return( f_end.lab ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( f_end.lab ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override { return( f_end.lab ); }

/*--------------------------------------------------------------------------*/

 /// extends Solver::int_par_type_S with the ThermalUnitDPSolver parameters
 enum int_par_type_TUDPS {
  intParMinN = intLastAlgPar , ///< min time horizon for the parallel DP path
  /**< The per-ON-node Economic Dispatch sweep in compute_EDPs() is run in
   * parallel (over the FastFlow workers set by intMaxThread) only when the
   * time horizon is at least intParMinN; below it the thread-dispatch
   * overhead is not amortised and the solver stays serial. Defaults to the
   * compile-time TUDPS_PAR_MIN_N; set it to 0 to force the parallel path on
   * every horizon. */
  intLastAlgParTUDPS ///< 1st allowed new int parameter for derived classes
  };

 using Solver::set_par;  // keep the other set_par() overloads visible

 /// honoured parameters: intMaxThread (workers) and intParMinN (threshold)
 void set_par( idx_type par , int value ) override {
  if( par == intMaxThread ) { f_max_thread = value; return; }
  if( par == intParMinN ) { f_par_min_n = value; return; }
  Solver::set_par( par , value );
  }

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == intMaxThread ) return( f_max_thread );
  if( par == intParMinN ) return( f_par_min_n );
  return( Solver::get_int_par( par ) );
  }

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( Solver::get_num_int_par() + intLastAlgParTUDPS - intLastAlgPar );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  return( par == intParMinN ? TUDPS_PAR_MIN_N
                            : Solver::get_dflt_int_par( par ) );
  }

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  return( name == "intParMinN" ? intParMinN
                               : Solver::int_par_str2idx( name ) );
  }

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string name = "intParMinN";
  return( idx == intParMinN ? name : Solver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 /// builds the graph
 void build_graph( void );

 /// computes the EDPs
 void compute_EDPs( void );

 /// implements the min-path algorithm
 void min_path( void );

 /// computes the variable values and the total cost
 void compute_solutions( void );

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

 /// stage of the computation
 enum stage_value
 {
  start = 0 ,
  graph_OK = 1 ,
  edps_OK = 2 ,
  path_OK = 3 ,
  sol_OK = 4
 };

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS EDSolver --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// base class for the Economic Dispatch Solver
 /** EDSolver is a base class that defines a minimal interface between the
  * ThermalUnitDPSolver and the solvers of the individual Economic Dispatch
  * Problems that give the cost of the arc in the DP. This is geared towards
  * solvers that can cheaply compute all the costs of all the arcs
  * ( h , h ), ( h , h + 1 ), ..., ( h , n ) in one blow, as the DP
  * solver does. However, it being virtual other implementations may be
  * considered. */

 class EDSolver
 {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

  EDSolver( Index h , ThermalUnitDPSolver * s )
   : f_h( h ) , f_solver( s ) {}

  virtual ~EDSolver() = default;

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

  /// compute the cost vector z_h[ h ], ..., z_h[ n - 1 ], z_h[ d ]
  /** Compute the costs of the arcs ( h , h ), ( h , h + 1 ), ...
   * ( h , n - 1 ), ( h , d ), where t is the end of the horizon and d is the
   * end node. The picture for h == 2 and n == 6 is
   *
   *       (0,1)   (1,1)   (2,1) -------+-------+-------+
   *                             \       \       \       \-> (t)
   *       (0,0)   (1,0)   (2,0)   (3,0)   (4,0)   (5,0)
   *
   * That is, these are t - h [ = 6 - 2 = 4 ] values corresponding to the
   * costs of the arcs in the DP graph of type ( h, h + 1 ), ( h, h + 2 ),
   * ..., ( h, n - 1 ), and finally the special arc ( h, d ) [ ( 2, 3 ),
   * ( 2, 4 ), ( 2, 5 ), ( 2, t ) ]. These are written in the positions h,
   * h + 1, ..., t - 1 [ 2 , 3 , 4 , 5 ] of the vector cost.
   *
   * The cost of each arc ( h , k ) for h < k <= n - 1 is ED( h , k - 1 ),
   * corresponding to the fact that the unit remains on from h to k - 1
   * included, but it is off at k. The cost of the special arc ( h , t )
   * corresponds to a "special" ED( h , n - 1 ) in which the unit remains
   * on from h to the end of the time horizon, comprised the last instant.
   * The difference is that in this last ED we do *not* assume the unit will
   * be shut down at n, as this is outside of the time horizon and whatever
   * happens to the unit then is of no concern here.
   *
   * More specifically, the point is that, due to the ramp-down constraints,
   * if the unit has to be down at time k, then it must enter in a "shutdown
   * trajectory" in the previous time instants, so that the power at k - 1
   * is the right one to stop. This constrains the ED, resulting in a higher
   * cost. This means that forcing the shut down at the end of n - 1 is never
   * economical: it is in principle better to allow the unit do what it wants
   * (which may comprise autonomously entering in a shutdown trajectory if
   * this is the optimal thing to do, as this is not prohibited). This is
   * why the constraints of the "special" ED( h , n - 1 ) do not include the
   * one forcing the power of the unit at the last time instant to be the
   * shutdown one, unlike for all the other ED( h , k ). */

  virtual void compute_costs( std::vector< double > & costs ) = 0;

/*--------------------------------------------------------------------------*/
  /// compute optimal power values p_h[ h ], p_h[ h + 1 ], ..., p_h[ k - 1 ]
  /** After compute_costs() have been called once, it is possible to call
   * compute_power_variables( k ) for h <= k <= t - 1 to get the optimal
   * power values corresponding to the arc ( h , k - 1 ); note that this also
   * works for the arc ( h , d ) by passing k = t, as we cheat so that the two
   * correspond to the same ED (see compute_costs()). The optimal values of
   * the power variables are written in the positions h, h + 1, ..., k - 1 of
   * the vector p. The solution depends on k, but this method is typically
   * only called for one particular value of k >= h during the final
   * computation of the optimal solution to the whole 1UC. */

  virtual void compute_power_variables( Index k ,
                                        std::vector< double > & p ) = 0;

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

  /// the ThermalUnitDPSolver using this EDSolver
  ThermalUnitDPSolver * f_solver;

  Index f_h;  ///< the initial time instant for this EDSolver

 };  // end( class( EDSolver ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DPEDSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// class solving the Economic Dispatch problem via Dynamic Programming
 /** DPEDSolver derives from EDSolver and solves the Economic Dispatch
  * problem by means of a Dynamic Programming approach. */

 class DPEDSolver : public EDSolver
 {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

  DPEDSolver( Index h , ThermalUnitDPSolver * s );

  virtual ~DPEDSolver() = default;

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

  void compute_costs( std::vector< double > & costs ) override;

  void compute_power_variables( Index k ,
                                std::vector< double > & p ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

  /// coefficients for a variable of the objective function
  struct coeff_t {
   double alfa;
   double beta;
   double gamma;
  };

  /// cost coefficients of the objective function
  std::vector< coeff_t > coeffs;

  /// indices for a piece of the (piece-wise) objective function
  struct pos_t {
   int begt;
   int begm;
  };

  /** For each k = h, ..., n - 1 the vector contains the indices of the pieces
   * of the objective function. */
  std::vector< pos_t > pos;

  /// unconstrained optimal power values
  std::vector< double > unc_p;

  /// constrained optimal power values
  std::vector< double > con_p;

  std::vector< double > m;
  std::vector< int > v;

  /// ping-pong half-size multiplier: 1 normally, larger when reserves are
  /// present so the augmented (g_t-split) pieces fit in each half
  Index f_rmul{ 1 };

  /// economic-dispatch sweep with the residual-ramp reserve reward folded in
  /** Alternative to compute_costs() taken when some reserve is rewarded.
   * Instead of the single-parabola coeffs[] / m[] sweep + capacity-band g
   * add, it carries a convex piecewise-quadratic value function
   * \f$ z_{h,k} \f$ and, at each interior step, replaces the pure-energy
   * ramp-window projection with the base class' sliding_min_corr(), the
   * deliverability-correct residual-ramp transition. The energy-only path
   * (compute_costs()) is untouched. Fills the same costs[ k ] = ED( h , k )
   * arc costs. */
  void compute_costs_reserve( std::vector< double > & costs );

 };  // end( class( DPEDSolver ) );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 class node;  // forward declaration of node

/*--------------------------------------------------------------------------*/
 /// an arc

 class arc
 {
  public:

  arc( void ) : cost1( 0 ) , cost2( 0 ) , tail( nullptr ) {}

  ~arc() = default;

  double cost1;  ///< the non-power-dependent part of the cost (fixed, SUC)
  double cost2;  ///< the power-dependent part of the cost
  node * tail;   ///< (pointer to) the tail node

 };  // end( class( arc ) )

/*--------------------------------------------------------------------------*/
 /// a node

 class node
 {
  public:

  node( void ) : lab( 0 ) , pred( nullptr ) , DPS( nullptr ) {}

  ~node() = default;

  double lab;                 ///< the label of the node
  node * pred;                ///< the predecessor of the node in the path
  EDSolver * DPS;             ///< the ED solver of the node (non-owning: the
                              ///< solvers are owned by the solver's pool)
  std::vector< arc > v_arcs;  ///< the Forward Star of the node

 };  // end( class( node ) )

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 Index h_of_node( node * n ) {
  if( n == &f_start )  // the source should be "-1", but we make it 0
   return( 0 );
  if( n == &f_end )    // the destination
   return( time_horizon );
  if( n->DPS )          // an ON-node
   return( n - v_on_nodes.data() );
  // else it must be an OFF-node, this is never called on the destination
  return( n - v_off_nodes.data() );
 }

/*--------------------------------------------------------------------------*/

 // reset label and predecessor of a node

 static void init_node( node & nde ) {
  nde.lab = TUDPINF;
  nde.pred = nullptr;
 }

/*--------------------------------------------------------------------------*/

 // reset a node for a fresh build_graph(), reusing its arc storage and the
 // pooled ED solvers; lab == 0 means "not yet proved reachable from s"

 static void reset_node( node & nde ) {
  nde.lab = 0;
  nde.pred = nullptr;
  nde.DPS = nullptr;
  nde.v_arcs.clear();  // keeps the allocated capacity
 }

/*--------------------------------------------------------------------------*/

 // return the pooled ED solver for the ON node at instant i, allocating it
 // the first time (and reusing it, buffers included, on later re-solves)

 DPEDSolver * get_on_ed( Index i );

/*--------------------------------------------------------------------------*/

 // do the scanning of the forward star of a node

 static void process_node( node & nde ) {
  for( auto a : nde.v_arcs ) {
   const auto nl = nde.lab + a.cost1 + a.cost2;
   if( ( a.tail )->lab > nl ) {
    ( a.tail )->lab = nl;
    ( a.tail )->pred = &nde;
   }
  }
 }

/*--------------------------------------------------------------------------*/

 void load_parameters( void );

/*--------------------------------------------------------------------------*/

 /// read the fixed status of the Variable of the ThermalUnitBlock
 /** Reads which commitment (and design) Variable are fixed, translating the
  * commitment fixings into the nxt_off / nxt_on tables that build_graph()
  * uses to prune the incompatible arcs; throws if any Variable that the DP
  * cannot honor (active power, reserves, start-up, shut-down, reactive) is
  * fixed. */
 void load_fixings( void );

/*--------------------------------------------------------------------------*/

 void process_modifications( void );

 // returns true if everything needs to be reset
 bool guts_of_process_modifications( const p_Mod mod );

/*--------------------------------------------------------------------------*/

 /// true iff some reserve price is negative, i.e. the reserve may be rewarded
 /** When true the economic dispatch takes the residual-ramp reserve path
  * (compute_costs_reserve() / the eff_disc* tables); when false it is the
  * plain quadratic sweep. The reserve model itself (reserve_alloc /
  * reserve_alloc_band / reserve_reward / build_reserve_discount /
  * sliding_min_corr / ...) lives in the base class ThermalUnitDPSolverBase
  * and is inherited. */
 void recover_schedule( std::vector< double > & p ,
                        std::vector< double > & u ,
                        std::vector< double > & pr ,
                        std::vector< double > & sr ,
                        std::vector< double > & q ,
                        bool & built ) const;
 ///< recovers the schedule the DP has found: power, commitment, reserves,
 ///< reactive power, and whether the unit is built

 bool reserve_rewarded( void ) const;

/*--------------------------------------------------------------------------*/

 double compute_startup_costs( Index h , Index k ) {
  // hook point for a time-dependent SUC( h , k ) formula
  if( startup_costs.empty() )
   return( 0 );
  return( startup_costs[ k ] );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 // NOTE: the data loaded from the ThermalUnitBlock (power/ramp/cost bounds,
 // the spinning-reserve factors and prices, the design and commitment-fixing
 // data, the scalar parameters), the PieceQuad / PQFun value-function type
 // and its operations, and the whole reserve model (reserve_alloc /
 // reserve_alloc_band / reserve_reward / build_reserve_discount /
 // reserve_corr_argmin / sliding_min_corr / ...) live in the base class
 // ThermalUnitDPSolverBase and are inherited. Only the solver-specific
 // reactive increment, the on/off graph, the ED-solver pool and the output
 // are declared here.

 // reactive power q[t] in [reactive_min[t] + reactive_min_on[t] u[t],
 // reactive_max[t] + reactive_max_on[t] u[t]], separable from the
 // active-power DP but gated by the commitment u[t]; reactive_linear_term[t]
 // is its dualized linear cost (empty if the unit has no reactive power).
 // The off-box reward is priced as the constant Q_star; when the box is
 // gated (the _on vectors are non-empty) the per-on-period increment
 // reactive_delta[t] = r_on - r_off is added to the fixed cost of every
 // on-period. fill_reactive_delta() (re)builds reactive_delta from the
 // current reactive price and boxes.
 std::vector< double > reactive_delta;

 void fill_reactive_delta( void );

 /// per-period reserve discount \f$ g_t( p ) \f$ as convex piecewise-linear
 /// PQFuns, precomputed once per compute_EDPs() (empty unless reserve is
 /// rewarded)
 /** eff_disc[ t ] is the interior variant (band cap max_power[ t ]);
  * eff_disc_su[ t ] the start-up variant (cap bound_on[ t ]), added at the
  * first period of an on-interval where there is no predecessor (hence no
  * ramp coupling). The shut-down variant is handled by recomputing the
  * closing transition with the shut-down cap inside compute_costs_reserve(),
  * so no eff_disc_sd table is needed. Built by the base class'
  * build_reserve_discount(). */
 std::vector< PQFun > eff_disc;
 std::vector< PQFun > eff_disc_su;

 char stage;                       ///< what has been computed

 node f_start;                     ///< starting node
 node f_end;                       ///< ending node

 std::vector< node > v_on_nodes;   ///< vector of ON nodes
 std::vector< node > v_off_nodes;  ///< vector of OFF nodes

 /// pool of ED solvers, one slot per time instant, reused across re-solves
 /// so that the (dominant) cost of allocating their O(n) buffers is paid
 /// only once per time horizon rather than at every structural re-solve
 std::vector< std::unique_ptr< DPEDSolver > > v_on_eds;
 std::unique_ptr< DPEDSolver > f_start_ed;  ///< ED solver of the source when
                                            ///< it acts as an ON node
 Index ed_pool_th{ 0 };            ///< time horizon the ED pool was built for

 int f_max_thread{ 0 };           ///< intMaxThread: 0 = all cores, 1 = serial
 int f_par_min_n{ TUDPS_PAR_MIN_N }; ///< intParMinN: min horizon for parallel

#if TUDPS_PARALLEL
 /// FastFlow parallel-for engine for compute_EDPs(), one per solver instance
 /// (created on first parallel use); its worker threads are reused across
 /// re-solves, and the per-worker cost scratch lives in f_tcost
 std::unique_ptr< ff::ParallelFor > f_pf;
 std::vector< std::vector< double > > f_tcost;
#endif

 std::vector< double > P;          ///< power values
 std::vector< bool > U;            ///< commitment values

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( ThermalUnitDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File ThermalUnitDPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
